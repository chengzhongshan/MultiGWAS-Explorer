#!/usr/bin/env perl
use strict;
use warnings;

use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use Getopt::Long qw(GetOptions);
use Text::CSV;

my ($candidates, $pfile, $bfile, $plink2, $populations, $output_leads,
    $output_audit, $output_cache, $output_status, $snp_column, $chr_column,
    $bp_column, $signal_column, $signal_columns);
my $signal_threshold = 1;
my $r2_threshold = 0.1;
my $window_kb = 1000;
my $threads = 4;
my $population_rule = 'ANY';
my $max_leads = 0;

$plink2 = $ENV{PLINK2} || 'plink2';
$populations = 'EUR EAS';
$snp_column = 'SNP';
$chr_column = 'CHR';
$bp_column = 'BP';
$signal_column = '';
$signal_columns = '';

GetOptions(
    'candidates=s'       => \$candidates,
    'pfile=s'            => \$pfile,
    'bfile=s'            => \$bfile,
    'plink2=s'           => \$plink2,
    'populations=s'      => \$populations,
    'population-rule=s'  => \$population_rule,
    'snp-column=s'       => \$snp_column,
    'chr-column=s'       => \$chr_column,
    'bp-column=s'        => \$bp_column,
    'signal-column=s'    => \$signal_column,
    'signal-columns=s'   => \$signal_columns,
    'signal-threshold=f' => \$signal_threshold,
    'r2-threshold=f'     => \$r2_threshold,
    'window-kb=f'        => \$window_kb,
    'threads=i'          => \$threads,
    'max-leads=i'        => \$max_leads,
    'output-leads=s'     => \$output_leads,
    'output-audit=s'     => \$output_audit,
    'output-cache=s'     => \$output_cache,
    'output-status=s'    => \$output_status,
) or die usage();

die usage() unless defined($candidates) && (defined($pfile) xor defined($bfile))
    && defined($output_leads) && defined($output_audit)
    && (length($signal_column) xor length($signal_columns));
die "--r2-threshold must be between 0 and 1\n"
    unless $r2_threshold >= 0 && $r2_threshold <= 1;
$population_rule = uc($population_rule);
die "--population-rule must be ANY or ALL\n"
    unless $population_rule =~ /\A(?:ANY|ALL)\z/;
my @pops = grep { length } map { uc($_) } split /[\s,]+/, $populations;
die "At least one population is required\n" unless @pops;

my ($rows, $header, $sep) = load_candidates();
die "No candidates passed the signal filter\n" unless @$rows;
my %by_snp = map { $_->{key} => $_ } @$rows;

my $tmp = tempdir('plink2_1kg_greedy_clump.XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $ids = File::Spec->catfile($tmp, 'candidates.ids.txt');
open my $id_fh, '>:raw', $ids or die "Cannot write $ids: $!\n";
print {$id_fh} "$_->{snp}\n" for @$rows;
close $id_fh;
my $subset = File::Spec->catfile($tmp, 'candidate_reference');
my @extract = plink_input();
push @extract, '--extract', native_path($ids), '--snps-only', 'just-acgt',
    '--max-alleles', 2, '--make-pgen', '--threads', $threads,
    '--out', native_path($subset);
run(@extract);

my %frequency;
for my $pop (@pops) {
    my $prefix = File::Spec->catfile($tmp, lc($pop) . '.freq');
    run(plink_exe(), '--pfile', native_path($subset), '--keep-if', 'SuperPop',
        '==', $pop, '--freq', 'counts', '--threads', $threads,
        '--out', native_path($prefix));
    $frequency{$pop} = read_acount("$prefix.acount");
}
write_reference_status(\%frequency);

my (%pruned, %audit, @selected);
open my $cache_fh, '>:raw', $output_cache
    or die "Cannot write $output_cache: $!\n" if defined $output_cache;
print {$cache_fh} join("\t", qw(query_snp proxy_snp ld_population proxy_r2 source)), "\n"
    if $cache_fh;

for my $lead (@$rows) {
    next if $pruned{$lead->{key}};
    last if $max_leads > 0 && @selected >= $max_leads;
    push @selected, $lead;
    my $lead_rank = scalar @selected;
    my (%proxy, @estimable_pops);
    for my $pop (@pops) {
        next unless polymorphic($frequency{$pop}{$lead->{snp}});
        push @estimable_pops, $pop;
        print {$cache_fh} join("\t", $lead->{key}, $lead->{key}, $pop, 1,
            'PLINK2_1KG_DIRECT'), "\n" if $cache_fh;
        my $prefix = File::Spec->catfile($tmp, sprintf('lead_%05d_%s', $lead_rank, lc($pop)));
        run(plink_exe(), '--pfile', native_path($subset), '--keep-if', 'SuperPop',
            '==', $pop, '--mac', 1, '--ld-snp', $lead->{snp},
            '--r2-phased', 'allow-ambiguous-allele', 'cols=chrom,pos,id',
            '--ld-window-kb', $window_kb, '--ld-window-r2', $r2_threshold,
            '--threads', $threads, '--out', native_path($prefix), '--silent');
        my $report = -s "$prefix.vcor" ? "$prefix.vcor" : "$prefix.vcor.zst";
        my $edges = read_vcor($report, $lead->{snp});
        for my $edge (@$edges) {
            my ($other, $r2) = @$edge;
            my $key = uc($other);
            next unless exists $by_snp{$key};
            next if $by_snp{$key}{rank} <= $lead->{rank} || $pruned{$key};
            $proxy{$key}{pops}{$pop} = 1;
            if (!defined($proxy{$key}{r2}) || $r2 > $proxy{$key}{r2}) {
                $proxy{$key}{r2} = $r2;
                $proxy{$key}{max_pop} = $pop;
            }
            print {$cache_fh} join("\t", $lead->{key}, $key, $pop, $r2,
                'PLINK2_1KG_DIRECT'), "\n" if $cache_fh;
        }
    }

    my $query_status = @estimable_pops == @pops ? 'ESTIMABLE'
                     : @estimable_pops ? 'PARTIAL_ESTIMABLE' : 'NOT_ESTIMABLE';
    $audit{$lead->{key}} = {
        action => 'SELECTED_LEAD', lead => $lead->{key}, lead_rank => $lead_rank,
        method => @estimable_pops ? 'LD' : 'UNRESOLVED_NO_LD',
        query_status => $query_status, population => '', r2 => '',
    };
    for my $key (keys %proxy) {
        next if $population_rule eq 'ALL' && keys(%{$proxy{$key}{pops}}) < @pops;
        $pruned{$key} = 1;
        my @edge_pops = sort keys %{$proxy{$key}{pops}};
        $audit{$key} = {
            action => 'PRUNED_LD', lead => $lead->{key}, lead_rank => $lead_rank,
            method => 'LD', query_status => '', r2 => $proxy{$key}{r2},
            population => @edge_pops > 1 ? join('+', @edge_pops) : $edge_pops[0],
        };
    }
}
close $cache_fh if $cache_fh;
for my $row (@$rows) {
    $audit{$row->{key}} ||= { action => 'NOT_REACHED', lead => '', lead_rank => '',
        method => '', query_status => '', population => '', r2 => '' };
}
write_leads($header, $sep, \@selected, \%audit);
write_audit($rows, \%audit);
print "TOP_HIT_LD_SOURCE\tPLINK2_1KG_DIRECT\n";
print "CANDIDATES\t", scalar(@$rows), "\nSELECTED_LEADS\t", scalar(@selected), "\n";
print "PRUNED_LD\t", scalar(grep { $_->{action} eq 'PRUNED_LD' } values %audit), "\n";
print "NOT_ESTIMABLE_LEADS\t", scalar(grep { $_->{action} eq 'SELECTED_LEAD' && $_->{query_status} eq 'NOT_ESTIMABLE' } values %audit), "\n";

sub load_candidates {
    open my $fh, '<:raw', $candidates or die "Cannot read $candidates: $!\n";
    my $first = <$fh> // die "Empty candidate file: $candidates\n";
    seek($fh, 0, 0);
    my $sep = index($first, "\t") >= 0 ? "\t" : ',';
    my $csv = Text::CSV->new({ binary => 1, auto_diag => 2, sep_char => $sep });
    my $header = $csv->getline($fh);
    my %idx = map { uc(trim($header->[$_])) => $_ } 0 .. $#$header;
    for my $name ($snp_column, $chr_column, $bp_column) {
        die "Candidate file lacks $name\n" unless exists $idx{uc $name};
    }
    my @signals = length($signal_column) ? ($signal_column) : split(/[\s,]+/, $signal_columns);
    for my $name (@signals) {
        die "Candidate file lacks $name\n" unless exists $idx{uc $name};
    }
    my (@rows, %seen);
    while (my $v = $csv->getline($fh)) {
        my $snp = trim($v->[$idx{uc $snp_column}]);
        next unless length $snp;
        my $key = uc($snp);
        next if $seen{$key}++;
        my @signals = grep { defined } map { numeric($v->[$idx{uc $_}]) } @signals;
        next unless @signals;
        my $signal = $signals[0];
        for my $candidate_signal (@signals) {
            $signal = $candidate_signal if $candidate_signal < $signal;
        }
        next unless $signal > 0 && $signal <= $signal_threshold;
        my $bp = numeric($v->[$idx{uc $bp_column}]);
        next unless defined $bp;
        push @rows, { values => $v, snp => $snp, key => $key,
            chr => trim($v->[$idx{uc $chr_column}]), bp => $bp, signal => $signal };
    }
    close $fh;
    @rows = sort { $a->{signal} <=> $b->{signal} || $a->{chr} cmp $b->{chr}
        || $a->{bp} <=> $b->{bp} || $a->{key} cmp $b->{key} } @rows;
    $rows[$_]{rank} = $_ + 1 for 0 .. $#rows;
    return (\@rows, $header, $sep);
}

sub read_acount {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "Cannot read $path: $!\n";
    my @h = split /\t/, scalar(<$fh>), -1;
    s/^#// for @h; s/[\r\n]+\z// for @h;
    my %i = map { uc($h[$_]) => $_ } 0 .. $#h;
    my %result;
    while (my $line = <$fh>) {
        $line =~ s/[\r\n]+\z//;
        my @f = split /\t/, $line, -1;
        $result{$f[$i{ID}]} = [0 + $f[$i{ALT_CTS}], 0 + $f[$i{OBS_CT}]];
    }
    close $fh;
    return \%result;
}

sub polymorphic {
    my ($counts) = @_;
    return $counts && $counts->[1] > 0 && $counts->[0] > 0 && $counts->[0] < $counts->[1];
}

sub read_vcor {
    my ($path, $lead) = @_;
    my $fh;
    if ($path =~ /\.zst\z/) { open $fh, '-|', ($ENV{ZSTD} || 'zstd'), '-dc', $path or die $!; }
    else { open $fh, '<:raw', $path or die "Cannot read $path: $!\n"; }
    my @h = split /\t/, scalar(<$fh>), -1;
    s/^#// for @h; s/[\r\n]+\z// for @h;
    my %i = map { uc($h[$_]) => $_ } 0 .. $#h;
    my $rcol = exists($i{PHASED_R2}) ? 'PHASED_R2' : 'R2';
    my @result;
    while (my $line = <$fh>) {
        $line =~ s/[\r\n]+\z//;
        my @f = split /\t/, $line, -1;
        my ($a, $b, $r2) = @f[$i{ID_A}, $i{ID_B}, $i{$rcol}];
        my $other = uc($a) eq uc($lead) ? $b : uc($b) eq uc($lead) ? $a : undef;
        push @result, [$other, 0 + $r2] if defined($other) && numeric($r2) >= $r2_threshold;
    }
    close $fh;
    return \@result;
}

sub write_reference_status {
    my ($frequency) = @_;
    return unless defined $output_status;
    open my $fh, '>:raw', $output_status or die "Cannot write $output_status: $!\n";
    print {$fh} join("\t", qw(variant ld_population reference_status alt_count obs_count maf source)), "\n";
    for my $row (@$rows) { for my $pop (@pops) {
        my $c = $frequency->{$pop}{$row->{snp}};
        my ($status, $alt, $obs, $maf) = ('ABSENT_REFERENCE', '', '', '');
        if ($c) { ($alt, $obs) = @$c; $status = polymorphic($c) ? 'POLYMORPHIC' : $obs ? 'MONOMORPHIC' : 'NO_GENOTYPE_OBSERVATIONS';
            if ($obs) { $maf = $alt / $obs; $maf = 1 - $maf if $maf > .5; } }
        print {$fh} join("\t", $row->{key}, $pop, $status, $alt, $obs,
            length($maf) ? sprintf('%.8g', $maf) : '', 'PLINK2_1KG_DIRECT'), "\n";
    }}
    close $fh;
}

sub write_leads {
    my ($header, $input_sep, $selected, $audit) = @_;
    my $sep = $output_leads =~ /\.tsv\z/i ? "\t" : ',';
    my $csv = Text::CSV->new({ binary => 1, eol => "\n", sep_char => $sep });
    open my $fh, '>:raw', $output_leads or die "Cannot write $output_leads: $!\n";
    $csv->print($fh, [@$header, qw(INDEPENDENCE_METHOD LD_SOURCE LD_POPULATIONS LD_POPULATION_RULE LD_QUERY_STATUS LD_LEAD_RANK LD_R2_THRESHOLD LD_WINDOW_KB)]);
    for my $row (@$selected) { my $a = $audit->{$row->{key}};
        $csv->print($fh, [@{$row->{values}}, $a->{method}, 'PLINK2_1KG_DIRECT', join(' ', @pops),
            $population_rule, $a->{query_status}, $a->{lead_rank}, $r2_threshold, $window_kb]); }
    close $fh;
}

sub write_audit {
    my ($rows, $audit) = @_;
    open my $fh, '>:raw', $output_audit or die "Cannot write $output_audit: $!\n";
    print {$fh} join("\t", qw(CHR BP candidate_snp lead_snp selection_action independence_method ld_source query_status ld_population candidate_rank signal lead_rank prune_r2)), "\n";
    for my $row (@$rows) { my $a = $audit->{$row->{key}};
        print {$fh} join("\t", $row->{chr}, $row->{bp}, $row->{key}, $a->{lead}, $a->{action},
            $a->{method}, 'PLINK2_1KG_DIRECT', $a->{query_status}, $a->{population}, $row->{rank},
            $row->{signal}, $a->{lead_rank}, $a->{r2}), "\n"; }
    close $fh;
}

sub plink_input {
    my @cmd = (plink_exe());
    if (defined $pfile) { push @cmd, '--pfile', native_path($pfile); push @cmd, 'vzs' if -e "$pfile.pvar.zst"; }
    else { push @cmd, '--bfile', native_path($bfile); }
    return @cmd;
}
sub plink_exe { return $plink2; }
sub run { my @cmd = @_; system(@cmd); die "Command failed (exit $?): @cmd\n" if $? != 0; }
sub numeric { my $v = trim($_[0]); return undef unless $v =~ /\A[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\z/; return 0 + $v; }
sub trim { my $v = defined($_[0]) ? $_[0] : ''; $v =~ s/^\s+|\s+$//g; return $v; }
sub native_path {
    my $path = shift; my $abs = abs_path($path) || File::Spec->rel2abs($path);
    if ($abs =~ m{^/mnt/([A-Za-z])/(.*)$}) { my ($d,$r)=(uc($1),$2); $r =~ s{/}{\\}g; return "$d:\\$r"; }
    if ($^O =~ /cygwin/i && $abs =~ m{^/}) { if (open my $fh, '-|', 'cygpath', '-w', $abs) { my $w=<$fh>//''; close $fh; $w =~ s/[\r\n]+\z//; return $w if length $w; } }
    return $abs;
}
sub usage { return <<'USAGE';
Usage: clump_top_hits_with_plink2_1kg.pl --candidates FILE
  (--pfile PREFIX | --bfile PREFIX) --output-leads FILE --output-audit FILE
  (--signal-column P | --signal-columns "P1 P2") [--output-cache FILE]
  [--output-status FILE] [--populations "EUR EAS"] [--population-rule ANY]
  [--r2-threshold 0.1] [--window-kb 1000] [--plink2 EXE]
USAGE
}
