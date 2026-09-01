#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw($GunzipError);
use File::Basename qw(dirname);
use File::Path qw(make_path);

sub usage {
    return <<"USAGE";
Usage:
  perl pdl_gunplot_manhattan.pl --data wide.tsv.gz --out-prefix prefix [options]

Options:
  --data FILE.tsv.gz        Input wide GWAS TSV.gz.
  --out-prefix PREFIX       Output prefix for .png/.gp/.plot.tsv/.manifest.tsv.
  --pcols A,B,C             Comma-separated P columns.
  --labels A|B|C            Pipe-separated track labels.
  --title TEXT              Optional plot title. Default: none.
  --width N                 PNG width. Default: 1800
  --height N                PNG height. Default: 1250
  --top-logp FLOAT          Per-track cap. Default: 12
  --sig FLOAT               Significance threshold. Default: 1e-6
  --min-logp FLOAT          Drop points below this. Default: 0.5
  --keep-all-logp FLOAT     Keep all points above this. Default: 3
  --thin-mod N              Background thinning factor. Default: 40
  --remove-x-chr            Exclude chromosome X / 23 from the rendered figure.
  --gnuplot PATH            gnuplot executable. Default: gnuplot
USAGE
}

my %opt = (
    width         => 1800,
    height        => 1250,
    top_logp      => 12,
    sig           => '1e-6',
    min_logp      => 0.5,
    keep_all_logp => 3,
    thin_mod      => 40,
    remove_x_chr  => 0,
    gnuplot       => 'gnuplot',
    title         => '',
);

GetOptions(
    'data=s'         => \$opt{data},
    'out-prefix=s'   => \$opt{out_prefix},
    'pcols=s'        => \$opt{pcols},
    'labels=s'       => \$opt{labels},
    'title=s'        => \$opt{title},
    'width=i'        => \$opt{width},
    'height=i'       => \$opt{height},
    'top-logp=f'     => \$opt{top_logp},
    'sig=s'          => \$opt{sig},
    'min-logp=f'     => \$opt{min_logp},
    'keep-all-logp=f'=> \$opt{keep_all_logp},
    'thin-mod=i'     => \$opt{thin_mod},
    'remove-x-chr!'  => \$opt{remove_x_chr},
    'gnuplot=s'      => \$opt{gnuplot},
) or die usage();

die usage() unless $opt{data} && $opt{out_prefix};
die "Input file not found: $opt{data}\n" unless -s $opt{data};

my @pcols = grep { length } map { trim($_) } split /,/, ($opt{pcols} // '');
die "--pcols is required\n" unless @pcols;

my @labels = map { trim($_) } split /\|/, ($opt{labels} // '');
@labels = map { nice_label($_) } @pcols unless @labels == @pcols;

my ($header, $idx_ref) = read_header($opt{data});
my %idx = %{$idx_ref};
my @resolved_pcols;
for my $need (@pcols) {
    my $resolved = resolve_existing_column($need, \%idx);
    die "Required column missing from input: $need\n" unless defined $resolved;
    push @resolved_pcols, $resolved;
}
for my $need ('CHR', 'BP') {
    die "Required column missing from input: $need\n" unless exists $idx{$need};
}

my $plot_tsv = "$opt{out_prefix}.plot.tsv";
my $gp_file  = "$opt{out_prefix}.gp";
my $png_file = "$opt{out_prefix}.png";
my $manifest = "$opt{out_prefix}.manifest.tsv";

my $chr_max = scan_chr_layout(
    data => $opt{data},
    idx  => \%idx,
    remove_x_chr => $opt{remove_x_chr},
);

my @chrs = sort { chr_order($a) <=> chr_order($b) } keys %{$chr_max};
die "No valid chromosome rows found in $opt{data}\n" unless @chrs;
my ($chr_offset, $chr_mid, $max_x) = chromosome_offsets($chr_max, \@chrs);

my ($rows_kept, $rows_scanned, $rows_thinned, $track_points) = write_plot_tsv(
    data       => $opt{data},
    idx        => \%idx,
    pcols      => \@resolved_pcols,
    plot_tsv   => $plot_tsv,
    top_logp   => $opt{top_logp},
    min_logp   => $opt{min_logp},
    keep_all   => $opt{keep_all_logp},
    thin_mod   => $opt{thin_mod},
    chr_offset => $chr_offset,
    remove_x_chr => $opt{remove_x_chr},
);

my $sig_y = safe_neglog10($opt{sig});
write_gnuplot(
    gp_file   => $gp_file,
    plot_tsv  => $plot_tsv,
    png_file  => $png_file,
    title     => $opt{title},
    width     => $opt{width},
    height    => $opt{height},
    top_logp  => $opt{top_logp},
    sig_y     => $sig_y,
    labels    => \@labels,
    chrs      => \@chrs,
    chr_mid   => $chr_mid,
    max_x     => $max_x,
);

system($opt{gnuplot}, $gp_file) == 0
    or die "gnuplot failed for $gp_file\n";

open my $mf, '>', $manifest or die "Cannot write $manifest: $!\n";
print {$mf} join("\t", qw(METRIC VALUE)), "\n";
print {$mf} join("\t", 'input', $opt{data}), "\n";
print {$mf} join("\t", 'png', $png_file), "\n";
print {$mf} join("\t", 'gnuplot_script', $gp_file), "\n";
print {$mf} join("\t", 'plot_tsv', $plot_tsv), "\n";
print {$mf} join("\t", 'title', $opt{title}), "\n";
print {$mf} join("\t", 'pcols', join(',', @pcols)), "\n";
print {$mf} join("\t", 'resolved_pcols', join(',', @resolved_pcols)), "\n";
print {$mf} join("\t", 'labels', join('|', @labels)), "\n";
print {$mf} join("\t", 'remove_x_chr', ($opt{remove_x_chr} ? 1 : 0)), "\n";
print {$mf} join("\t", 'rows_scanned', $rows_scanned), "\n";
print {$mf} join("\t", 'rows_kept', $rows_kept), "\n";
print {$mf} join("\t", 'rows_thinned', $rows_thinned), "\n";
for my $i (0 .. $#pcols) {
    print {$mf} join("\t", "points_$pcols[$i]", ($track_points->{$i} || 0)), "\n";
}
close $mf or die "Cannot close $manifest: $!\n";

print "PNG\t$png_file\n";
print "MANIFEST\t$manifest\n";

sub read_header {
    my ($path) = @_;
    my $fh = IO::Uncompress::Gunzip->new($path)
        or die "Cannot read $path: $GunzipError\n";
    my $line = <$fh>;
    close $fh;
    die "Input is empty: $path\n" unless defined $line;
    chomp $line;
    $line =~ s/\r$//;
    my @cols = split /\t/, $line, -1;
    my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
    return ($line, \%idx);
}

sub scan_chr_layout {
    my (%args) = @_;
    my $path      = $args{data};
    my $idx       = $args{idx};
    my $remove_x_chr = $args{remove_x_chr} ? 1 : 0;
    my $chr_idx = $idx->{CHR};
    my $bp_idx  = $idx->{BP};

    my %chr_max;
    my $fh = IO::Uncompress::Gunzip->new($path)
        or die "Cannot read $path: $GunzipError\n";
    <$fh>;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        my $chr = normalize_chr($f[$chr_idx]);
        my $bp  = numeric($f[$bp_idx]);
        next if $remove_x_chr && is_chr_x($chr);
        next unless length($chr) && defined $bp && $bp > 0;
        $chr_max{$chr} = $bp if !exists $chr_max{$chr} || $bp > $chr_max{$chr};
    }
    close $fh;
    return \%chr_max;
}

sub write_plot_tsv {
    my (%args) = @_;
    my $path       = $args{data};
    my $idx        = $args{idx};
    my $pcols      = $args{pcols};
    my $plot_tsv   = $args{plot_tsv};
    my $top_logp   = $args{top_logp};
    my $min_logp   = $args{min_logp};
    my $keep_all   = $args{keep_all};
    my $thin_mod   = $args{thin_mod};
    my $chr_offset = $args{chr_offset};
    my $remove_x_chr = $args{remove_x_chr} ? 1 : 0;

    my $chr_idx = $idx->{CHR};
    my $bp_idx  = $idx->{BP};
    my @track_points = (0) x scalar(@{$pcols});
    my $rows_scanned = 0;
    my $rows_thinned = 0;
open my $out, '>', $plot_tsv or die "Cannot write $plot_tsv: $!\n";
print {$out} join("\t", qw(X Y TRACK CHR BP LOGP CHR_COLOR_INDEX)), "\n";
    my $fh = IO::Uncompress::Gunzip->new($path)
        or die "Cannot read $path: $GunzipError\n";
    <$fh>;
    my $rows_kept = 0;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        $rows_scanned++;
        my @f = split /\t/, $line, -1;
        my $chr = normalize_chr($f[$chr_idx]);
        my $bp  = numeric($f[$bp_idx]);
        next if $remove_x_chr && is_chr_x($chr);
        next unless length($chr) && defined $bp && $bp > 0 && exists $chr_offset->{$chr};
        for my $track_i (0 .. $#{$pcols}) {
            my $p_col = $pcols->[$track_i];
            my $p = numeric($f[ $idx->{$p_col} ]);
            next unless defined $p && $p > 0 && $p <= 1;
            my $logp = safe_neglog10($p);
            next unless defined $logp && $logp >= $min_logp;
            if ($logp < $keep_all && $thin_mod > 1) {
                my $hash = simple_hash(join(':', $chr, $bp, $track_i));
                if ($hash % $thin_mod != 0) {
                    $rows_thinned++;
                    next;
                }
            }
            my $x = $chr_offset->{$chr} + $bp;
            my $capped = $logp > $top_logp ? $top_logp : $logp;
            my $y = $track_i * $top_logp + $capped;
            my $chr_color_idx = chr_palette_index($chr);
            print {$out} join(
                "\t",
                sprintf('%.4f', $x),
                sprintf('%.4f', $y),
                $track_i,
                $chr,
                $bp,
                sprintf('%.6f', $logp),
                $chr_color_idx,
            ), "\n";
            $track_points[$track_i]++;
            $rows_kept++;
        }
    }
    close $fh;
    close $out or die "Cannot close $plot_tsv: $!\n";

    my %track_count = map { $_ => $track_points[$_] } 0 .. $#track_points;
    return ($rows_kept, $rows_scanned, $rows_thinned, \%track_count);
}

sub chromosome_offsets {
    my ($chr_max, $chrs) = @_;
    my (%offset, %mid);
    my $cursor = 0;
    for my $chr (@{$chrs}) {
        $offset{$chr} = $cursor;
        $mid{$chr} = $cursor + ($chr_max->{$chr} || 1) / 2;
        $cursor += ($chr_max->{$chr} || 1) + 1_000_000;
    }
    return (\%offset, \%mid, $cursor);
}

sub write_gnuplot {
    my (%args) = @_;
    open my $gp, '>', $args{gp_file} or die "Cannot write $args{gp_file}: $!\n";
    my @xtics = map { sprintf('"%s" %.0f', $_, $args{chr_mid}{$_}) } @{ $args{chrs} };
    my @ytics = repeated_panel_ytics_manhattan(scalar(@{ $args{labels} }), $args{top_logp});

    print {$gp} "set terminal png size $args{width},$args{height}\n";
    print {$gp} "set output '" . escape_gp($args{png_file}) . "'\n";
    print {$gp} "set datafile separator '\\t'\n";
    if (defined $args{title} && length $args{title}) {
        print {$gp} "set title \"" . escape_gp($args{title}) . "\"\n";
    }
    else {
        print {$gp} "unset title\n";
    }
    print {$gp} "set xlabel 'Chromosome'\n";
    print {$gp} "set ylabel '-log10(P)'\n";
    print {$gp} "unset key\n";
    print {$gp} "set border 3\n";
    print {$gp} "set lmargin 10\n";
    print {$gp} "set rmargin 3\n";
    print {$gp} "set tics out nomirror\n";
    print {$gp} "set grid ytics lc rgb '#dddddd' dt 2\n";
    print {$gp} "set xrange [0:$args{max_x}]\n";
    print {$gp} "set yrange [0:" . ((scalar(@{ $args{labels} }) * $args{top_logp}) + 1.2) . "]\n";
    print {$gp} "set xtics (" . join(', ', @xtics) . ")\n";
    print {$gp} "set ytics (" . join(', ', @ytics) . ")\n";
    print {$gp} "set mytics 2\n";
    print {$gp} "set palette maxcolors 24 defined (";
    my @colors = sas_chr_palette();
    print {$gp} join(', ', map { sprintf('%d "%s"', $_ + 1, $colors[$_]) } 0 .. $#colors) . ")\n";
    print {$gp} "unset colorbox\n";

    my $arrow_id = 1;
    for my $i (0 .. $#{ $args{labels} }) {
        my $base = $i * $args{top_logp};
        print {$gp} "set arrow $arrow_id from 0," . ($base + $args{sig_y}) . " to $args{max_x}," . ($base + $args{sig_y}) . " nohead dt 2 lc rgb '#777777' lw 1\n";
        $arrow_id++;
        next unless $i > 0;
        print {$gp} "set arrow $arrow_id from 0,$base to $args{max_x},$base nohead lc rgb '#bbbbbb' lw 1\n";
        $arrow_id++;
    }
    for my $i (0 .. $#{ $args{labels} }) {
        my $panel_y = $i * $args{top_logp} + $args{top_logp} - 0.95;
        my $xmid = $args{max_x} / 2;
        print {$gp} "set label " . (100 + $i) . " \"" . escape_gp($args{labels}[$i]) . "\" at $xmid,$panel_y center font ',18'\n";
    }

    print {$gp} "plot '" . escape_gp($args{plot_tsv}) . "' using 1:2:7 every ::1 with points pt 7 ps 0.33 lc palette\n";
    close $gp or die "Cannot close $args{gp_file}: $!\n";
}

sub repeated_panel_ytics_manhattan {
    my ($track_count, $top_logp) = @_;
    my @ticks;
    for my $i (0 .. ($track_count - 1)) {
        my $base = $i * $top_logp;
        for my $tick (0, 2, 4, 6, 8) {
            next if $tick > $top_logp;
            push @ticks, sprintf('"%d" %.3f', $tick, $base + $tick);
        }
    }
    return @ticks;
}

sub numeric {
    my ($x) = @_;
    return undef unless defined $x && $x ne '';
    return undef if $x =~ /^(?:NA|NaN|null|\.)$/i;
    return undef unless $x =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
    return 0 + $x;
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

sub safe_neglog10 {
    my ($value) = @_;
    my $num = numeric($value);
    return undef unless defined $num && $num > 0;
    return -log($num) / log(10);
}

sub simple_hash {
    my ($text) = @_;
    my $hash = 2166136261;
    for my $c (unpack('C*', $text)) {
        $hash ^= $c;
        $hash = ($hash * 16777619) & 0xffffffff;
    }
    return $hash;
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
    if ($requested =~ /^(.*)_STD_Z$/) {
        push @cand, "${1}_STD_DIFF_Z", "${1}_DIFF_Z";
    }
    if ($requested =~ /^(.*)_FEMALE_P$/) {
        push @cand, "${1}_GROUP1_P";
    }
    if ($requested =~ /^(.*)_MALE_P$/) {
        push @cand, "${1}_GROUP2_P";
    }
    if ($requested =~ /^(.*)_FEMALE_Z$/) {
        push @cand, "${1}_GROUP1_Z";
    }
    if ($requested =~ /^(.*)_MALE_Z$/) {
        push @cand, "${1}_GROUP2_Z";
    }
    return @cand;
}

sub nice_label {
    my ($text) = @_;
    $text =~ s/_/ /g;
    $text =~ s/\bSTD\b/standardized/ig;
    $text =~ s/\bDIFF\b/diff/ig;
    $text =~ s/\bP\b/P/ig;
    $text =~ s/\s+/ /g;
    $text =~ s/^\s+|\s+$//g;
    return $text;
}

sub trim {
    my ($x) = @_;
    $x //= '';
    $x =~ s/^\s+//;
    $x =~ s/\s+$//;
    return $x;
}

sub escape_gp {
    my ($x) = @_;
    $x //= '';
    $x =~ s{\\}{\\\\}g;
    $x =~ s{'}{''}g;
    $x =~ s{"}{\\"}g;
    return $x;
}

sub sas_chr_palette {
    return (
        '#0072bd', '#d95319', '#edb120', '#7e2f8e',
        '#77ac30', '#4dbeee', '#a2142f',
    );
}

sub chr_palette_index {
    my ($chr) = @_;
    my $ord = chr_order($chr);
    my @palette = sas_chr_palette();
    $ord = 1 unless defined $ord && $ord > 0 && $ord < 10_000;
    return (($ord - 1) % scalar(@palette)) + 1;
}
