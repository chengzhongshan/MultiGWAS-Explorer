#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);

my $input =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.tsv.gz';
my $out_prefix =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.manhattan';
my $group_col       = 'PAIR_TAG';
my $p_cols          = 'STD_DIFF_P,DIFF_P';
my $chr_col         = 'CHR';
my $pos_col         = 'BP';
my $snp_col         = 'SNP';
my $width           = 1800;
my $height          = 1200;
my $top_logp        = 10;
my $gwas_threshold  = 7.30103;
my $min_logp        = 0.5;
my $keep_all_logp   = 3;
my $thin_mod        = 10;
my $chr_gap         = 5_000_000;
my $point_size      = 0.22;
my $gnuplot         = 'gnuplot';
my $title           = 'PGC schizophrenia sex differential GWAS Manhattan Plot';
my $reuse_plot_tsv  = 0;
my $remove_x_chr    = 0;

GetOptions(
    'input=s'          => \$input,
    'out-prefix=s'     => \$out_prefix,
    'group-col=s'      => \$group_col,
    'p-cols=s'         => \$p_cols,
    'chr-col=s'        => \$chr_col,
    'pos-col=s'        => \$pos_col,
    'snp-col=s'        => \$snp_col,
    'width=i'          => \$width,
    'height=i'         => \$height,
    'top-logp=f'       => \$top_logp,
    'gwas-threshold=f' => \$gwas_threshold,
    'min-logp=f'       => \$min_logp,
    'keep-all-logp=f'  => \$keep_all_logp,
    'thin-mod=i'       => \$thin_mod,
    'chr-gap=i'        => \$chr_gap,
    'point-size=f'     => \$point_size,
    'gnuplot=s'        => \$gnuplot,
    'title=s'          => \$title,
    'reuse-plot-tsv!'  => \$reuse_plot_tsv,
    'remove-x-chr!'    => \$remove_x_chr,
) or die usage();

die "Input file not found: $input\n" unless -s $input;
die "--thin-mod must be at least 1\n" unless $thin_mod >= 1;

my @p_cols = grep { $_ ne '' } map { trim($_) } split /,/, $p_cols;
die "No P-value columns supplied with --p-cols\n" unless @p_cols;

my $plot_tsv = "$out_prefix.plot.tsv";
my $gp_file  = "$out_prefix.gnuplot";
my $png_file = "$out_prefix.png";
my $manifest = "$out_prefix.manifest.tsv";

my ($header, $idx) = read_header($input);
for my $required ($chr_col, $pos_col, $group_col, @p_cols) {
    die "Required column $required not found in input header\n" unless exists $idx->{$required};
}
my $has_snp_col = exists $idx->{$snp_col};

my ($chr_max, $groups) = scan_layout($input, $idx);
my @chrs = sort { chr_order($a) <=> chr_order($b) } keys %$chr_max;
my @groups = sort {
    canonical_group_rank($a) <=> canonical_group_rank($b)
      || uc($a) cmp uc($b)
} keys %$groups;
die "No chromosome positions found in $input\n" unless @chrs;
die "No groups found in $input column $group_col\n" unless @groups;

my (@tracks, %track_i);
for my $p_col (@p_cols) {
    for my $group (@groups) {
        push @tracks, { p_col => $p_col, group => $group, label => short_track_label($group, $p_col) };
    }
}
@tracks = sort {
    track_sort_rank($a->{group}, $a->{p_col}) <=> track_sort_rank($b->{group}, $b->{p_col})
      || $a->{label} cmp $b->{label}
} @tracks;
for my $i (0 .. $#tracks) {
    $track_i{ $tracks[$i]{p_col} }{ $tracks[$i]{group} } = $i;
}

my ($chr_offset, $chr_mid, $max_x) = chromosome_offsets($chr_max, \@chrs);
my $plot_counts = ($reuse_plot_tsv && -s $plot_tsv)
    ? rewrite_plot_data_from_existing($plot_tsv, \%track_i)
    : write_plot_data($input, $idx, \@p_cols, \%track_i, $chr_offset, $plot_tsv);
write_gnuplot($gp_file, $plot_tsv, $png_file, $chr_mid, \@chrs, \@tracks, $max_x);

system($gnuplot, $gp_file) == 0 or die "gnuplot failed for $gp_file\n";
write_manifest($manifest, $plot_tsv, $gp_file, $png_file, \@tracks, $plot_counts, $chr_max, $max_x);

print "Input:    $input\n";
print "Plot TSV: $plot_tsv\n";
print "Gnuplot:  $gp_file\n";
print "PNG:      $png_file\n";
print "Manifest: $manifest\n";
print "Tracks:   ", scalar(@tracks), "\n";
print "Points:   $plot_counts->{kept}\n";

sub read_header {
    my ($path) = @_;
    open my $fh, '-|', "zcat '$path'" or die "Cannot read $path with zcat: $!\n";
    my $h = <$fh>;
    close $fh;
    die "Input is empty: $path\n" unless defined $h;
    chomp $h;
    $h =~ s/\r$//;
    $h =~ s/^#//;
    my @cols = split /\t/, $h, -1;
    my %idx;
    for my $i (0 .. $#cols) {
        $idx{$cols[$i]} = $i;
    }
    return ($h, \%idx);
}

sub scan_layout {
    my ($path, $idx) = @_;
    open my $fh, '-|', "zcat '$path'" or die "Cannot read $path with zcat: $!\n";
    <$fh>;
    my (%chr_max, %groups);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        my @v = split /\t/, $line, -1;
        my $chr = normalize_chr($v[ $idx->{$chr_col} ]);
        my $bp  = numeric($v[ $idx->{$pos_col} ]);
        my $grp = $v[ $idx->{$group_col} ] // '';
        next unless defined $bp && $bp > 0 && $chr ne '' && $grp ne '';
        next if $remove_x_chr && is_chr_x($chr);
        $chr_max{$chr} = $bp if !exists $chr_max{$chr} || $bp > $chr_max{$chr};
        $groups{$grp} = 1;
    }
    close $fh;
    return (\%chr_max, \%groups);
}

sub rewrite_plot_data_from_existing {
    my ($path, $track_i) = @_;
    my $tmp = "$path.restyled";
    open my $in, '<', $path or die "Cannot read $path: $!\n";
    open my $out, '>', $tmp or die "Cannot write $tmp: $!\n";
    my $header = <$in>;
    die "Existing plot TSV is empty: $path\n" unless defined $header;
    print {$out} $header;
    chomp $header;
    $header =~ s/\r$//;
    my @cols = split /\t/, $header, -1;
    my %col = map { $cols[$_] => $_ } 0 .. $#cols;
    for my $required (qw(X LOGP P GROUP P_COL SNP CHR_NUM)) {
        die "Existing plot TSV lacks column $required\n" unless exists $col{$required};
    }
    my %counts = (rows => 0, candidates => 0, kept => 0, skipped => 0, reused => 1);
    while (my $line = <$in>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        $counts{rows}++;
        my @f = split /\t/, $line, -1;
        my $group = $f[ $col{GROUP} ];
        my $p_col = $f[ $col{P_COL} ];
        my $track = $track_i->{$p_col}{$group};
        next unless defined $track;
        my $logp = numeric($f[ $col{LOGP} ]);
        next unless defined $logp;
        my $capped = $logp > $top_logp ? $top_logp : $logp;
        my $y = $track * $top_logp + $capped;
        print {$out} join(
            "\t",
            $f[ $col{X} ],
            $y,
            $f[ $col{CHR_NUM} ],
            $track,
            $f[ $col{LOGP} ],
            $f[ $col{P} ],
            $group,
            $p_col,
            $f[ $col{SNP} ],
        ), "\n";
        $counts{kept}++;
    }
    close $in;
    close $out;
    rename $tmp, $path or die "Cannot replace $path with restyled plot TSV: $!\n";
    return \%counts;
}

sub count_plot_rows {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    my $count = -1;
    $count++ while <$fh>;
    close $fh;
    return $count < 0 ? 0 : $count;
}

sub chromosome_offsets {
    my ($chr_max, $chrs) = @_;
    my (%offset, %mid);
    my $cursor = 0;
    for my $chr (@$chrs) {
        $offset{$chr} = $cursor;
        $mid{$chr} = $cursor + $chr_max->{$chr} / 2;
        $cursor += $chr_max->{$chr} + $chr_gap;
    }
    return (\%offset, \%mid, $cursor);
}

sub write_plot_data {
    my ($path, $idx, $p_cols, $track_i, $chr_offset, $outpath) = @_;
    open my $fh, '-|', "zcat '$path'" or die "Cannot read $path with zcat: $!\n";
    open my $out, '>', $outpath or die "Cannot write $outpath: $!\n";
    print {$out} join("\t", qw(X Y CHR_NUM TRACK_INDEX LOGP P GROUP P_COL SNP)), "\n";
    <$fh>;

    my %counts = (rows => 0, candidates => 0, kept => 0, skipped => 0);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        my @v = split /\t/, $line, -1;
        $counts{rows}++;

        my $chr = normalize_chr($v[ $idx->{$chr_col} ]);
        my $bp  = numeric($v[ $idx->{$pos_col} ]);
        my $grp = $v[ $idx->{$group_col} ] // '';
        next if $remove_x_chr && is_chr_x($chr);
        next unless defined $bp && exists $chr_offset->{$chr} && exists $track_i->{$p_cols->[0]}{$grp};
        my $snp = $has_snp_col ? ($v[ $idx->{$snp_col} ] // '') : '';
        my $x = $chr_offset->{$chr} + $bp;
        my $chr_num = chr_order($chr);

        for my $p_col (@$p_cols) {
            my $p = numeric($v[ $idx->{$p_col} ]);
            next unless defined $p && $p > 0 && $p <= 1;
            my $logp = -log($p) / log(10);
            next if $logp < $min_logp;
            $counts{candidates}++;
            if ($logp < $keep_all_logp && $thin_mod > 1) {
                my $h = simple_hash(join(':', $chr, $bp, $grp, $p_col, $snp));
                if ($h % $thin_mod != 0) {
                    $counts{skipped}++;
                    next;
                }
            }
            my $track = $track_i->{$p_col}{$grp};
            my $capped = $logp > $top_logp ? $top_logp : $logp;
            my $y = $track * $top_logp + $capped;
            print {$out} join("\t", $x, $y, $chr_num, $track, fmt($logp), $p, $grp, $p_col, $snp), "\n";
            $counts{kept}++;
        }
    }
    close $fh;
    close $out;
    return \%counts;
}

sub write_gnuplot {
    my ($gp, $data, $png, $chr_mid, $chrs, $tracks, $max_x) = @_;
    open my $out, '>', $gp or die "Cannot write $gp: $!\n";

    my @xtics = map { sprintf('"%s" %.0f', $_, $chr_mid->{$_}) } @$chrs;
    my @ytics = repeated_panel_ytics($tracks);
    my $label_x = $max_x / 2;

    print {$out} "set terminal png size $width,$height\n";
    print {$out} "set output '$png'\n";
    print {$out} "set datafile separator '\\t'\n";
    print {$out} "set title '" . escape_gp_label($title) . "'\n";
    print {$out} "set xlabel 'Chromosome'\n";
    print {$out} "set ylabel '-Log10(p)'\n";
    print {$out} "unset key\n";
    print {$out} "set border 3\n";
    print {$out} "set lmargin 10\n";
    print {$out} "set rmargin 3\n";
    print {$out} "set tics out nomirror\n";
    print {$out} "set xrange [0:$max_x]\n";
    print {$out} "set yrange [0:", scalar(@$tracks) * $top_logp, "]\n";
    print {$out} "set xtics (", join(', ', @xtics), ")\n";
    print {$out} "set ytics (", join(', ', @ytics), ")\n";
    print {$out} "set mytics 2\n";
    print {$out} "set palette maxcolors 24 defined (";
    my @colors = sas_chr_palette();
    my @palette;
    for my $i (0 .. $#colors) {
        push @palette, sprintf('%d "%s"', $i + 1, $colors[$i]);
    }
    print {$out} join(', ', @palette), ")\n";
    print {$out} "unset colorbox\n";
    print {$out} "set grid ytics lc rgb '#dddddd'\n";

    my $arrow = 1;
    for my $i (0 .. $#$tracks) {
        my $base = $i * $top_logp;
        my $thr_y = $base + $gwas_threshold;
        print {$out} "set arrow $arrow from 0,$thr_y to $max_x,$thr_y nohead dt 2 lc rgb '#777777' lw 1\n";
        $arrow++;
        if ($i > 0) {
            print {$out} "set arrow $arrow from 0,$base to $max_x,$base nohead lc rgb '#bdbdbd' lw 1\n";
            $arrow++;
        }
        my $panel_y = $base + $top_logp - 1.2;
        print {$out} "set label ", ($i + 1), " \"" . escape_gp_label($tracks->[$i]{label}) . "\" at $label_x,$panel_y center font ',20'\n";
    }

    print {$out} "plot '$data' using 1:2:3 every ::1 with points pt 7 ps $point_size lc palette\n";
    close $out;
}

sub write_manifest {
    my ($path, $plot_tsv, $gp, $png, $tracks, $counts, $chr_max, $max_x) = @_;
    open my $out, '>', $path or die "Cannot write $path: $!\n";
    print {$out} join("\t", qw(METRIC VALUE)), "\n";
    print {$out} join("\t", 'input', $input), "\n";
    print {$out} join("\t", 'plot_tsv', $plot_tsv), "\n";
    print {$out} join("\t", 'gnuplot_script', $gp), "\n";
    print {$out} join("\t", 'png', $png), "\n";
    print {$out} join("\t", 'tracks', scalar(@$tracks)), "\n";
    print {$out} join("\t", 'rows_scanned', $counts->{rows}), "\n";
    print {$out} join("\t", 'candidate_points', $counts->{candidates}), "\n";
    print {$out} join("\t", 'points_plotted', $counts->{kept}), "\n";
    print {$out} join("\t", 'points_thinned', $counts->{skipped}), "\n";
    print {$out} join("\t", 'plot_tsv_reused', ($counts->{reused} ? 1 : 0)), "\n";
    print {$out} join("\t", 'min_logp', $min_logp), "\n";
    print {$out} join("\t", 'keep_all_logp', $keep_all_logp), "\n";
    print {$out} join("\t", 'thin_mod', $thin_mod), "\n";
    print {$out} join("\t", 'top_logp', $top_logp), "\n";
    print {$out} join("\t", 'gwas_threshold', $gwas_threshold), "\n";
    print {$out} join("\t", 'max_x', $max_x), "\n";
    print {$out} join("\t", 'track_labels', join(';', map { $_->{label} } @$tracks)), "\n";
    print {$out} "CHROMOSOME\tMAX_BP\n";
    for my $chr (sort { chr_order($a) <=> chr_order($b) } keys %$chr_max) {
        print {$out} join("\t", $chr, $chr_max->{$chr}), "\n";
    }
    close $out;
}

sub normalize_chr {
    my ($chr) = @_;
    return '' unless defined $chr;
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

sub numeric {
    my ($x) = @_;
    return undef unless defined $x;
    return undef if $x eq '' || $x =~ /^(?:NA|NaN|null|\.)$/i;
    return undef unless $x =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
    return 0 + $x;
}

sub simple_hash {
    my ($s) = @_;
    my $h = 2166136261;
    for my $c (unpack('C*', $s)) {
        $h ^= $c;
        $h = ($h * 16777619) & 0xffffffff;
    }
    return $h;
}

sub trim {
    my ($x) = @_;
    $x =~ s/^\s+//;
    $x =~ s/\s+$//;
    return $x;
}

sub escape_gp_label {
    my ($x) = @_;
    $x =~ s/\\/\\\\/g;
    $x =~ s/"/\\"/g;
    return $x;
}

sub short_track_label {
    my ($group, $p_col) = @_;
    my $root = canonical_group_root($group);
    if ($p_col eq 'GROUP2_P') {
        return "$root MALE association P";
    }
    if ($p_col eq 'GROUP1_P') {
        return "$root FEMALE association P";
    }
    if ($p_col eq 'STD_DIFF_P') {
        my $display = $root eq 'ALL' ? 'All' : $root;
        return "$display standardized diff P";
    }
    if ($p_col eq 'DIFF_P') {
        my $display = $root eq 'ALL' ? 'All' : $root;
        return "$display diff P";
    }
    return "$root $p_col";
}

sub canonical_group_root {
    my ($group) = @_;
    my $g = uc($group // '');
    return 'EUR' if $g =~ /EUR/;
    return 'ASN' if $g =~ /ASN/;
    return 'ALL' if $g =~ /ALL/;
    return $g;
}

sub canonical_group_rank {
    my ($group) = @_;
    my $root = canonical_group_root($group);
    return 0 if $root eq 'EUR';
    return 1 if $root eq 'ASN';
    return 2 if $root eq 'ALL';
    return 9;
}

sub track_sort_rank {
    my ($group, $p_col) = @_;
    my $root = canonical_group_root($group);
    my %assoc_rank = ( ALL => 0, ASN => 1, EUR => 2 );
    my %diff_rank  = ( ALL => 0, EUR => 1, ASN => 2 );
    if ($p_col eq 'STD_DIFF_P') {
        return ($diff_rank{$root} // 9);
    }
    if ($p_col eq 'DIFF_P') {
        return 20 + ($diff_rank{$root} // 9);
    }
    if ($p_col eq 'GROUP1_P') {
        return 40 + ($assoc_rank{$root} // 9) * 2;
    }
    if ($p_col eq 'GROUP2_P') {
        return 40 + ($assoc_rank{$root} // 9) * 2 + 1;
    }
    return 999;
}

sub repeated_panel_ytics {
    my ($tracks) = @_;
    my @ticks;
    for my $i (0 .. $#$tracks) {
        my $base = $i * $top_logp;
        for my $tick (0, 2, 4, 6, 8) {
            next if $tick > $top_logp;
            push @ticks, sprintf('"%d" %.3f', $tick, $base + $tick);
        }
    }
    return @ticks;
}

sub sas_chr_palette {
    return (
        '#0072bd', '#d95319', '#edb120', '#7e2f8e',
        '#77ac30', '#4dbeee', '#a2142f',
    );
}

sub fmt {
    my ($x) = @_;
    return sprintf('%.10g', $x);
}

sub usage {
    return <<"USAGE";
Usage:
  perl gnuplot_multiple_manhattan.pl [options]

Default input:
  PGC_SCZ_female_vs_male_diff_effects.stdized.tsv.gz

Options:
  --input FILE.tsv.gz       Input differential GWAS table
  --out-prefix PREFIX       Prefix for .plot.tsv, .gnuplot, .png, .manifest.tsv
  --group-col NAME          Track grouping column. Default: PAIR_TAG
  --p-cols A,B              P columns to plot. Default: STD_DIFF_P,DIFF_P
  --chr-col NAME            Chromosome column. Default: CHR
  --pos-col NAME            Position column. Default: BP
  --snp-col NAME            SNP ID column. Default: SNP
  --width N                 PNG width. Default: 1800
  --height N                PNG height. Default: 1200
  --top-logp FLOAT          Cap each track's -log10(P). Default: 10
  --gwas-threshold FLOAT    Reference line. Default: 7.30103
  --min-logp FLOAT          Drop points below this -log10(P). Default: 0.5
  --keep-all-logp FLOAT     Keep every point at or above this -log10(P). Default: 3
  --thin-mod N              Keep 1/N lower-signal background points. Default: 10
  --point-size FLOAT        Gnuplot point size. Default: 0.22
  --title TEXT              Plot title
  --reuse-plot-tsv          Reuse existing PREFIX.plot.tsv and only restyle/redraw
  --remove-x-chr            Exclude chromosome X / 23 from the rendered plot
  --gnuplot PATH            gnuplot executable. Default: gnuplot
USAGE
}
