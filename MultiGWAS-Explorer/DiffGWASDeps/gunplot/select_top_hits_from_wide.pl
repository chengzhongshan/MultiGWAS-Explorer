#!/usr/bin/env perl
use strict;
use warnings;

use FindBin qw($Bin);
use lib "$Bin/..";
use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw($GunzipError);
use TopHitMAF qw(
  numeric
  format_num
  derive_effect_af
  maf_from_effect_af
  parse_population_map
  infer_population_codes_for_text
  load_gnomad_lookup
  lookup_gnomad_maf
);

sub usage {
    return <<"USAGE";
Usage:
  perl select_top_hits_from_wide.pl --input wide.tsv.gz --output hits.tsv [options]

Options:
  --mode differential|targeted   Default: differential
  --focus-pvar NAME              Required unless --target-snps is used.
  --focus-prefix NAME            Optional explicit pair prefix for MAF lookup.
  --thresholds A,B,C             Threshold ladder. Default: 1e-6
  --top-hit-dist-bp N            Distance pruning span. Default: 1e6
  --max-hits N                   Optional max selected loci.
  --target-snps A,B,C            Explicit target SNP order.
  --remove-x-chr                 Exclude chromosome X / 23 from selected hits.
  --maf-threshold NUM            Minimum allowed MAF for top-hit selection.
                                 Default: 0.01. Use 0 to disable the filter.
  --gnomad-freq-file FILE        Optional lightweight gnomAD lookup TSV/TSV.GZ.
  --gnomad-pop-map MAP           Optional token-to-gnomAD-pop map such as
                                 EUR=NFE,ASN=EAS,AFR=AFR,ALL=AF
USAGE
}

my %opt = (
    mode            => 'differential',
    thresholds      => '1e-6',
    top_hit_dist_bp => '1e6',
    max_hits        => 0,
    target_snps     => '',
    remove_x_chr    => 0,
    maf_threshold   => 0.01,
    gnomad_freq_file => '',
    gnomad_pop_map   => '',
    focus_prefix     => '',
);

GetOptions(
    'input=s'            => \$opt{input},
    'output=s'           => \$opt{output},
    'mode=s'             => \$opt{mode},
    'focus-pvar=s'       => \$opt{focus_pvar},
    'focus-prefix=s'     => \$opt{focus_prefix},
    'thresholds=s'       => \$opt{thresholds},
    'top-hit-dist-bp=s'  => \$opt{top_hit_dist_bp},
    'max-hits=i'         => \$opt{max_hits},
    'target-snps=s'      => \$opt{target_snps},
    'remove-x-chr!'      => \$opt{remove_x_chr},
    'maf-threshold=f'    => \$opt{maf_threshold},
    'gnomad-freq-file=s' => \$opt{gnomad_freq_file},
    'gnomad-pop-map=s'   => \$opt{gnomad_pop_map},
) or die usage();

die usage() unless $opt{input} && $opt{output};
die "Input file not found: $opt{input}\n" unless -s $opt{input};
die "--focus-pvar is required for non-targeted selection\n"
    if !$opt{target_snps} && !$opt{focus_pvar};

my @target_snps = grep { length } map { trim($_) } split /,/, ($opt{target_snps} // '');
my @thresholds = grep { defined && length } map { trim($_) } split /[,\s]+/, $opt{thresholds};
@thresholds = ('1e-6') unless @thresholds;

my $pop_map = parse_population_map($opt{gnomad_pop_map});
my $gnomad_lookup = load_gnomad_lookup(file => $opt{gnomad_freq_file});

my $fh = IO::Uncompress::Gunzip->new($opt{input})
    or die "Cannot read $opt{input}: $GunzipError\n";
my $header_line = <$fh>;
die "Input file is empty: $opt{input}\n" unless defined $header_line;
chomp $header_line;
$header_line =~ s/\r$//;
my @header = split /\t/, $header_line, -1;
my %idx = map { $header[$_] => $_ } 0 .. $#header;
for my $need (qw(CHR BP SNP)) {
    die "Required column missing from input: $need\n" unless exists $idx{$need};
}
my $focus_pvar = $opt{focus_pvar};
if (!$opt{target_snps}) {
    $focus_pvar = resolve_existing_column($opt{focus_pvar}, \%idx);
    die "Focus P column missing from input: $opt{focus_pvar}\n" unless defined $focus_pvar;
}
my $focus_prefix = $opt{focus_prefix} || derive_focus_prefix($focus_pvar);

my @rows;
my $maf_filtered = 0;
my $maf_unknown = 0;
while (my $line = <$fh>) {
    chomp $line;
    $line =~ s/\r$//;
    next unless length $line;
    my @f = split /\t/, $line, -1;
    my $chr = normalize_chr($f[ $idx{CHR} ]);
    my $bp  = numeric($f[ $idx{BP} ]);
    my $snp = $f[ $idx{SNP} ] // '';
    next unless length($chr) && defined $bp && length($snp);
    next if $opt{remove_x_chr} && is_chr_x($chr);

    my $row = {
        CHR => $chr,
        BP  => $bp,
        SNP => $snp,
        raw => \@f,
    };
    $row->{focus_signal} = numeric($f[ $idx{$focus_pvar} ]) if !$opt{target_snps};
    $row->{maf_info} = annotate_differential_maf(
        row           => $row,
        header        => \@header,
        idx           => \%idx,
        focus_prefix  => $focus_prefix,
        maf_threshold => $opt{maf_threshold},
        gnomad_lookup => $gnomad_lookup,
        pop_map       => $pop_map,
    );
    if (!@target_snps && $opt{maf_threshold} > 0) {
        if (!$row->{maf_info}{pass}) {
            $maf_filtered++;
            $maf_unknown++ if ($row->{maf_info}{maf_source} || '') eq 'UNKNOWN';
            next;
        }
    }
    push @rows, $row;
}
close $fh;

my @selected;
if (@target_snps) {
    my %by_snp = map { uc($_->{SNP}) => $_ } @rows;
    my $order = 0;
    for my $snp (@target_snps) {
        my $hit = $by_snp{ uc($snp) } or next;
        $order++;
        $hit->{hit_order} = $order;
        $hit->{focus_signal} = defined $hit->{focus_signal} ? $hit->{focus_signal} : '';
        push @selected, $hit;
    }
}
else {
    my $chosen_threshold = $thresholds[-1];
    for my $thr (@thresholds) {
        my @cand = grep {
            defined $_->{focus_signal} && $_->{focus_signal} > 0 && $_->{focus_signal} < (0 + $thr)
        } @rows;
        if (@cand) {
            $chosen_threshold = $thr;
            last;
        }
    }

    my @candidates = grep {
        defined $_->{focus_signal} && $_->{focus_signal} > 0 && $_->{focus_signal} < (0 + $chosen_threshold)
    } @rows;
    @candidates = sort {
        chr_order($a->{CHR}) <=> chr_order($b->{CHR})
            ||
        $a->{focus_signal} <=> $b->{focus_signal}
            ||
        $a->{BP} <=> $b->{BP}
    } @candidates;

    my $half_window = (0 + $opt{top_hit_dist_bp}) / 2;
    my %selected_by_chr;
    my $order = 0;
    for my $hit (@candidates) {
        my $keep = 1;
        for my $sel (@{ $selected_by_chr{ $hit->{CHR} } || [] }) {
            if ($hit->{BP} >= $sel->{BP} - $half_window && $hit->{BP} <= $sel->{BP} + $half_window) {
                $keep = 0;
                last;
            }
        }
        next unless $keep;
        $order++;
        $hit->{hit_order} = $order;
        push @selected, $hit;
        push @{ $selected_by_chr{ $hit->{CHR} } }, $hit;
        last if $opt{max_hits} && @selected >= $opt{max_hits};
    }
}

open my $out, '>', $opt{output} or die "Cannot write $opt{output}: $!\n";
print {$out} join("\t", qw(
  hit_order CHR BP SNP focus_signal
  selected_maf maf_source gwas_group1_maf gwas_group2_maf gwas_pair_maf_min
  gnomad_maf gnomad_pops maf_filter_decision maf_filter_reason
)), "\n";
for my $hit (@selected) {
    my $maf = $hit->{maf_info} || {};
    print {$out} join("\t",
        $hit->{hit_order} || '',
        $hit->{CHR},
        $hit->{BP},
        $hit->{SNP},
        defined $hit->{focus_signal} ? $hit->{focus_signal} : '',
        map { defined $_ ? $_ : '' } (
            format_num($maf->{selected_maf}),
            $maf->{maf_source},
            format_num($maf->{gwas_group1_maf}),
            format_num($maf->{gwas_group2_maf}),
            format_num($maf->{gwas_pair_maf_min}),
            format_num($maf->{gnomad_maf}),
            $maf->{gnomad_pops},
            $maf->{decision},
            $maf->{reason},
        ),
    ), "\n";
}
close $out or die "Cannot close $opt{output}: $!\n";

print "OUTPUT\t$opt{output}\n";
print "HITS\t" . scalar(@selected) . "\n";
print "MAF_FILTERED\t$maf_filtered\n";
print "MAF_UNKNOWN\t$maf_unknown\n";

sub annotate_differential_maf {
    my (%args) = @_;
    my $row = $args{row} || {};
    my $raw = $row->{raw} || [];
    my $idx = $args{idx} || {};
    my $prefix = $args{focus_prefix} || '';
    my $thr = defined $args{maf_threshold} ? $args{maf_threshold} : 0;

    my $g1_maf = compute_pair_group_maf($raw, $idx, $prefix, 1);
    my $g2_maf = compute_pair_group_maf($raw, $idx, $prefix, 2);
    my @gwas_mafs = grep { defined $_ } ($g1_maf, $g2_maf);
    my $gwas_pair_maf_min;
    if (@gwas_mafs) {
        @gwas_mafs = sort { $a <=> $b } @gwas_mafs;
        $gwas_pair_maf_min = $gwas_mafs[0];
        return {
            selected_maf       => $gwas_pair_maf_min,
            maf_source         => 'GWAS',
            gwas_group1_maf    => $g1_maf,
            gwas_group2_maf    => $g2_maf,
            gwas_pair_maf_min  => $gwas_pair_maf_min,
            decision           => ($thr > 0 && $gwas_pair_maf_min <= $thr) ? 'FILTERED' : 'PASS',
            reason             => ($thr > 0 && $gwas_pair_maf_min <= $thr)
                ? sprintf('GWAS pair minimum MAF %.6g <= %.6g', $gwas_pair_maf_min, $thr)
                : sprintf('GWAS pair minimum MAF %.6g > %.6g', $gwas_pair_maf_min, $thr),
            pass               => ($thr > 0 ? ($gwas_pair_maf_min > $thr ? 1 : 0) : 1),
        };
    }

    my @pop_codes = infer_population_codes_for_text($prefix, $args{pop_map});
    my $gnomad = lookup_gnomad_maf(
        lookup    => $args{gnomad_lookup},
        record    => $row,
        pop_codes => \@pop_codes,
    );
    if ($gnomad && defined $gnomad->{maf}) {
        my @pairs = map { $_ . '=' . format_num($gnomad->{pop_mafs}{$_}) } sort keys %{ $gnomad->{pop_mafs} || {} };
        return {
            selected_maf => $gnomad->{maf},
            maf_source   => 'GNOMAD',
            gnomad_maf   => $gnomad->{maf},
            gnomad_pops  => join('|', @pairs),
            decision     => ($thr > 0 && $gnomad->{maf} <= $thr) ? 'FILTERED' : 'PASS',
            reason       => ($thr > 0 && $gnomad->{maf} <= $thr)
                ? sprintf('gnomAD fallback MAF %.6g <= %.6g', $gnomad->{maf}, $thr)
                : sprintf('gnomAD fallback MAF %.6g > %.6g', $gnomad->{maf}, $thr),
            pass         => ($thr > 0 ? ($gnomad->{maf} > $thr ? 1 : 0) : 1),
        };
    }

    return {
        maf_source => 'UNKNOWN',
        decision   => ($thr > 0 ? 'FILTERED' : 'PASS'),
        reason     => 'No GWAS or gnomAD MAF was available for this candidate',
        pass       => ($thr > 0 ? 0 : 1),
    };
}

sub compute_pair_group_maf {
    my ($raw, $idx, $prefix, $group_num) = @_;
    return undef unless defined $prefix && length $prefix;
    my $fa = numeric(field($raw, $idx, "${prefix}_GROUP${group_num}_FRQ_A"));
    my $fu = numeric(field($raw, $idx, "${prefix}_GROUP${group_num}_FRQ_U"));
    my $eaf = derive_effect_af(frq_a => $fa, frq_u => $fu);
    return maf_from_effect_af($eaf);
}

sub field {
    my ($raw, $idx, $name) = @_;
    return '' unless exists $idx->{$name};
    return $raw->[ $idx->{$name} ] // '';
}

sub derive_focus_prefix {
    my ($focus_pvar) = @_;
    return '' unless defined $focus_pvar && length $focus_pvar;
    my $prefix = $focus_pvar;
    $prefix =~ s/_(?:STD_DIFF_P|STD_P|DIFF_P|GROUP1_P|GROUP2_P|P)$//;
    return $prefix;
}

sub resolve_existing_column {
    my ($requested, $idx) = @_;
    return $requested if exists $idx->{$requested};
    my @candidates = alias_candidates($requested);
    for my $cand (@candidates) {
        return $cand if exists $idx->{$cand};
    }
    return undef;
}

sub alias_candidates {
    my ($requested) = @_;
    my @cand;
    if ($requested =~ /^(.*)_STD_P$/) {
        push @cand, "${1}_STD_DIFF_P", "${1}_DIFF_P";
    }
    if ($requested =~ /^(.*)_FEMALE_P$/) {
        push @cand, "${1}_GROUP1_P";
    }
    if ($requested =~ /^(.*)_MALE_P$/) {
        push @cand, "${1}_GROUP2_P";
    }
    return @cand;
}

sub normalize_chr {
    my ($chr) = @_;
    return '' unless defined $chr;
    $chr =~ s/^\s+|\s+$//g;
    $chr =~ s/^chr//i;
    return 'X' if $chr =~ /^(?:23|X)$/i;
    return $chr if $chr =~ /^\d+$/;
    return '';
}

sub chr_order {
    my ($chr) = @_;
    return 23 if defined $chr && uc($chr) eq 'X';
    return $chr if defined $chr && $chr =~ /^\d+$/;
    return 10_000;
}

sub is_chr_x {
    my ($chr) = @_;
    return 0 unless defined $chr;
    return $chr =~ /^(?:X|23)$/i ? 1 : 0;
}

sub trim {
    my ($x) = @_;
    $x //= '';
    $x =~ s/^\s+//;
    $x =~ s/\s+$//;
    return $x;
}
