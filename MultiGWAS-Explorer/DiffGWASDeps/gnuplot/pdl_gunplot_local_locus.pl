#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw($GunzipError);
use File::Path qw(make_path);
use File::Basename qw(dirname basename);

sub usage {
    return <<"USAGE";
Usage:
  perl pdl_gunplot_local_locus.pl --data wide.tsv.gz --snp rsID --out-prefix prefix [options]

Options:
  --window-bp N            Required.
  --pcols A,B,C            Required.
  --zcols A,B,C            Optional z-score columns for color mapping.
  --labels A|B|C           Optional display labels.
  --label-snps A,B,C       Ordered SNPs to mark and label at the top. The first
                           SNP remains the locus center supplied by --snp.
  --ld-snps A,B,C          Optional SNPs in LD with the query SNP(s). Matching
                           scatter points are overlaid with the selected marker.
  --ld-marker-symbol NAME  star, plus, cross, circle, square, triangle, diamond
                           (default: star).
  --ld-marker-color COLOR  Named or #RRGGBB color (default: black).
  --ld-display-mode MODE   none, markers, heatmap, or both (default: none).
  --ld-r2-values MAP       Comma-separated SNP:r2 values for LD proxies.
  --ld-population POP      Population shown in the LD inset (default: EUR).
  --ld-heatmap-colors LIST Low-to-high #RRGGBB colors for the separate LD scale.
  --title TEXT             Optional title.
  --gtf FILE.tsv           Optional extracted GTF subset TSV.
  --width N                Default: 1500
  --height N               Default: 1000
  --top-logp FLOAT         Default: 8
  --sig FLOAT              Default: 1e-6
  --gnuplot PATH           Default: gnuplot
  --hide-y-axis            Hide repeated y-axis for non-left combined panels.
  --bottom-snp-label TEXT  Override bottom vertical SNP label.
  --bottom-gene-label TEXT Optional second bottom vertical gene label.
USAGE
}

my %opt = (
    width    => 1500,
    height   => 1000,
    top_logp => 8,
    sig      => '1e-6',
    gnuplot  => 'gnuplot',
    title    => '',
    hide_y_axis => 0,
    ld_marker_symbol => 'star',
    ld_marker_color  => 'black',
    ld_display_mode  => 'none',
    ld_population    => 'EUR',
    ld_heatmap_colors=> '#f7fbff,#6baed6,#54278f',
);

GetOptions(
    'data=s'       => \$opt{data},
    'snp=s'        => \$opt{snp},
    'out-prefix=s' => \$opt{out_prefix},
    'window-bp=s'  => \$opt{window_bp},
    'pcols=s'      => \$opt{pcols},
    'zcols=s'      => \$opt{zcols},
    'labels=s'     => \$opt{labels},
    'label-snps=s' => \$opt{label_snps},
    'ld-snps=s'    => \$opt{ld_snps},
    'ld-marker-symbol=s' => \$opt{ld_marker_symbol},
    'ld-marker-color=s'  => \$opt{ld_marker_color},
    'ld-display-mode=s'  => \$opt{ld_display_mode},
    'ld-r2-values=s'     => \$opt{ld_r2_values},
    'ld-population=s'    => \$opt{ld_population},
    'ld-heatmap-colors=s'=> \$opt{ld_heatmap_colors},
    'title=s'      => \$opt{title},
    'gtf=s'        => \$opt{gtf},
    'width=i'      => \$opt{width},
    'height=i'     => \$opt{height},
    'top-logp=f'   => \$opt{top_logp},
    'sig=s'        => \$opt{sig},
    'gnuplot=s'    => \$opt{gnuplot},
    'hide-y-axis!' => \$opt{hide_y_axis},
    'bottom-snp-label=s' => \$opt{bottom_snp_label},
    'bottom-gene-label=s' => \$opt{bottom_gene_label},
) or die usage();

die usage() unless $opt{data} && $opt{snp} && $opt{out_prefix} && $opt{window_bp} && $opt{pcols};
die "Input file not found: $opt{data}\n" unless -s $opt{data};
my $ld_point_type = ld_marker_point_type($opt{ld_marker_symbol});
die "--ld-marker-color must be a named color or #RRGGBB\n"
    unless $opt{ld_marker_color} =~ /^(?:#[0-9A-Fa-f]{6}|[A-Za-z][A-Za-z0-9_-]*)$/;
$opt{ld_display_mode} = lc(trim($opt{ld_display_mode} || 'none'));
die "--ld-display-mode must be none, markers, heatmap, or both\n"
    unless $opt{ld_display_mode} =~ /^(?:none|markers|heatmap|both)$/;
$opt{ld_population} = uc(trim($opt{ld_population} || 'EUR'));
die "--ld-population must be AFR, AMR, ASN, or EUR\n"
    unless $opt{ld_population} =~ /^(?:AFR|AMR|ASN|EUR)$/;
my @ld_heatmap_colors = map { lc(trim($_)) } split /,/, $opt{ld_heatmap_colors};
die "--ld-heatmap-colors requires at least two comma-separated #RRGGBB colors\n"
    unless @ld_heatmap_colors >= 2 && !grep { $_ !~ /^#[0-9a-f]{6}$/ } @ld_heatmap_colors;

my @pcols = grep { length } map { trim($_) } split /,/, $opt{pcols};
my @zcols = grep { length } map { trim($_) } split /,/, ($opt{zcols} // '');
my @labels = map { nice_label(trim($_)) } split /\|/, ($opt{labels} // '');
@labels = map { nice_label($_) } @pcols unless @labels == @pcols;
my @label_snps;
my %seen_label_snp;
for my $snp (split /,/, ($opt{label_snps} // $opt{snp})) {
    $snp = trim($snp);
    next unless length $snp;
    next if $seen_label_snp{lc $snp}++;
    push @label_snps, $snp;
}
@label_snps = ($opt{snp}) unless @label_snps;
unshift @label_snps, $opt{snp}
    unless grep { lc($_) eq lc($opt{snp}) } @label_snps;
my %is_label_snp = map { lc($_) => 1 } @label_snps;
my @ld_snps;
my %is_ld_snp;
my %ld_r2_for;
for my $snp (split /,/, ($opt{ld_snps} // '')) {
    $snp = trim($snp);
    next unless length $snp;
    next if $is_label_snp{lc $snp};
    next if $is_ld_snp{lc $snp}++;
    push @ld_snps, $snp;
}
for my $item (split /,/, ($opt{ld_r2_values} // '')) {
    my ($snp, $r2) = split /:/, $item, 2;
    $snp = trim($snp // '');
    next unless length($snp) && defined($r2)
        && $r2 =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i;
    $r2 = cap_num(0 + $r2, 0, 1);
    $ld_r2_for{lc $snp} = $r2
        if !exists($ld_r2_for{lc $snp}) || $r2 > $ld_r2_for{lc $snp};
}
my $has_gtf = $opt{gtf} && -s $opt{gtf};

my $plot_tsv = "$opt{out_prefix}.plot.tsv";
my $gene_tsv = "$opt{out_prefix}.genes.tsv";
my $gp_file  = "$opt{out_prefix}.gp";
my $png_file = "$opt{out_prefix}.png";
my $manifest = "$opt{out_prefix}.manifest.tsv";

my $fh = IO::Uncompress::Gunzip->new($opt{data})
    or die "Cannot read $opt{data}: $GunzipError\n";
my $header_line = <$fh>;
die "Input file is empty: $opt{data}\n" unless defined $header_line;
chomp $header_line;
$header_line =~ s/\r$//;
my @header = split /\t/, $header_line, -1;
my %idx = map { $header[$_] => $_ } 0 .. $#header;
my @resolved_pcols;
for my $need (@pcols) {
    my $resolved = resolve_existing_column($need, \%idx);
    die "Required column missing from input: $need\n" unless defined $resolved;
    push @resolved_pcols, $resolved;
}
my @resolved_zcols;
for my $need (@zcols) {
    my $resolved = resolve_existing_column($need, \%idx);
    push @resolved_zcols, (defined $resolved ? $resolved : $need);
}
for my $need ('CHR', 'BP', 'SNP') {
    die "Required column missing from input: $need\n" unless exists $idx{$need};
}

my ($target_chr, $target_bp, @window_rows);
my %requested_snp_position;
while (my $line = <$fh>) {
    chomp $line;
    $line =~ s/\r$//;
    next unless length $line;
    my @f = split /\t/, $line, -1;
    my $chr = normalize_chr($f[ $idx{CHR} ]);
    my $bp  = numeric($f[ $idx{BP} ]);
    my $snp = $f[ $idx{SNP} ] // '';
    next unless length($chr) && defined $bp && length($snp);
    if (!defined $target_chr && lc($snp) eq lc($opt{snp})) {
        $target_chr = $chr;
        $target_bp  = $bp;
    }
    if ($is_label_snp{lc $snp} && !exists $requested_snp_position{lc $snp}) {
        $requested_snp_position{lc $snp} = { snp => $snp, chr => $chr, bp => $bp };
    }
    push @window_rows, \@f;
}
close $fh;

die "Target SNP $opt{snp} was not found in $opt{data}\n" unless defined $target_chr && defined $target_bp;
my $window_bp = 0 + $opt{window_bp};
my $start = $target_bp - $window_bp;
$start = 1 if $start < 1;
my $end = $target_bp + $window_bp;

my @locus = grep {
    normalize_chr($_->[ $idx{CHR} ]) eq $target_chr
        &&
    defined numeric($_->[ $idx{BP} ])
        &&
    numeric($_->[ $idx{BP} ]) >= $start
        &&
    numeric($_->[ $idx{BP} ]) <= $end
} @window_rows;
die "No rows found in locus window for $opt{snp}\n" unless @locus;
my @target_markers;
for my $requested (@label_snps) {
    my $marker = $requested_snp_position{lc $requested};
    die "Label SNP $requested was not found in $opt{data}\n" unless $marker;
    die "Label SNP $requested is outside the plotted chromosome/window centered on $opt{snp}\n"
        unless $marker->{chr} eq $target_chr && $marker->{bp} >= $start && $marker->{bp} <= $end;
    push @target_markers, { snp => $requested, bp => $marker->{bp} };
}

open my $pt, '>', $plot_tsv or die "Cannot write $plot_tsv: $!\n";
print {$pt} join("\t", qw(BP Y TRACK LOGP IS_TARGET SNP COLORVAL IS_LD LD_R2 LD_RGB)), "\n";
my $kept_points = 0;
my $ld_points = 0;
my %found_ld_snp;
my $has_zcols = @resolved_zcols == @resolved_pcols ? 1 : 0;
for my $row (@locus) {
    my $bp  = numeric($row->[ $idx{BP} ]);
    my $snp = $row->[ $idx{SNP} ] // '';
    for my $track_i (0 .. $#resolved_pcols) {
        my $p = numeric($row->[ $idx{ $resolved_pcols[$track_i] } ]);
        next unless defined $p && $p > 0 && $p <= 1;
        my $logp = safe_neglog10($p);
        next unless defined $logp;
        my $capped = $logp > $opt{top_logp} ? $opt{top_logp} : $logp;
        my $y = $track_i * $opt{top_logp} + $capped;
        my $is_target = $is_label_snp{lc $snp} ? 1 : 0;
        my $is_ld = $is_ld_snp{lc $snp} ? 1 : 0;
        my $ld_r2 = $is_ld && exists($ld_r2_for{lc $snp}) ? $ld_r2_for{lc $snp} : -1;
        my $ld_rgb = $ld_r2 >= 0 ? ld_rgb_integer($ld_r2, \@ld_heatmap_colors) : 0;
        if ($is_ld) {
            $ld_points++;
            $found_ld_snp{lc $snp} = $snp;
        }
        my $colorval = $track_i + 1;
        if ($has_zcols) {
            my $z = extract_requested_numeric($resolved_zcols[$track_i], $row, \%idx);
            $colorval = defined $z ? cap_num($z, -8, 8) : 0;
        }
        print {$pt} join("\t", $bp, sprintf('%.4f', $y), $track_i, sprintf('%.6f', $logp), $is_target, $snp, sprintf('%.4f', $colorval), $is_ld, sprintf('%.4f', $ld_r2), $ld_rgb), "\n";
        $kept_points++;
    }
}
close $pt or die "Cannot close $plot_tsv: $!\n";

my $gene_rows = 0;
if ($has_gtf) {
    open my $gf, '<', $opt{gtf} or die "Cannot read $opt{gtf}: $!\n";
    my $gh = <$gf>;
    die "Empty GTF subset: $opt{gtf}\n" unless defined $gh;
    chomp $gh;
    my @gcols = split /\t/, $gh, -1;
    my %gidx = map { $gcols[$_] => $_ } 0 .. $#gcols;
    for my $need (qw(chr genesymbol gene st en type)) {
        die "Required GTF subset column missing: $need\n" unless exists $gidx{$need};
    }
    my (%genes, %gene_exons);
    open my $gt, '>', $gene_tsv or die "Cannot write $gene_tsv: $!\n";
    print {$gt} join("\t", qw(TYPE GENE START END LANE)), "\n";
    while (my $line = <$gf>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        my $type = lc(($f[ $gidx{type} ] // ''));
        next unless $type eq 'gene' || $type eq 'exon';
        my $chr = normalize_chr($f[ $gidx{chr} ]);
        next unless length($chr) && $chr eq $target_chr;
        my $gs = numeric($f[ $gidx{st} ]);
        my $ge = numeric($f[ $gidx{en} ]);
        next unless defined $gs && defined $ge;
        next if $ge < $start || $gs > $end;
        my $gene = $f[ $gidx{genesymbol} ] || $f[ $gidx{gene} ] || 'GENE';
        if ($type eq 'gene') {
            if (!exists $genes{$gene}) {
                $genes{$gene} = { start => $gs, end => $ge };
            }
            else {
                $genes{$gene}{start} = $gs if $gs < $genes{$gene}{start};
                $genes{$gene}{end} = $ge if $ge > $genes{$gene}{end};
            }
        }
        else {
            my $exon_key = join(':', $gs, $ge);
            $gene_exons{$gene}{$exon_key} = [$gs, $ge];
        }
    }
    my @lanes;
    my %gene_lane;
    my $plot_bp_span = max_num(1, $end - $start + 1);
    for my $gene (sort { $genes{$a}{start} <=> $genes{$b}{start} || $genes{$a}{end} <=> $genes{$b}{end} || $a cmp $b } keys %genes) {
        my ($lane_start, $lane_end) = lane_reservation_range(
            gene       => $gene,
            gene_start => $genes{$gene}{start},
            gene_end   => $genes{$gene}{end},
            plot_start => $start,
            plot_end   => $end,
            plot_width => ($opt{width} || 1500),
            plot_bp_span => $plot_bp_span,
        );
        my $lane = allocate_lane(\@lanes, $lane_start, $lane_end);
        $gene_lane{$gene} = $lane;
        print {$gt} join("\t", 'gene', $gene, $genes{$gene}{start}, $genes{$gene}{end}, $lane), "\n";
        for my $exon_key (sort {
            $gene_exons{$gene}{$a}[0] <=> $gene_exons{$gene}{$b}[0]
              || $gene_exons{$gene}{$a}[1] <=> $gene_exons{$gene}{$b}[1]
        } keys %{ $gene_exons{$gene} || {} }) {
            my ($es, $ee) = @{ $gene_exons{$gene}{$exon_key} };
            print {$gt} join("\t", 'exon', $gene, $es, $ee, $lane), "\n";
        }
        $gene_rows++;
    }
    close $gt or die "Cannot close $gene_tsv: $!\n";
    close $gf;
}

my $sig_y = safe_neglog10($opt{sig});
my $gene_height = $has_gtf ? max_num(6.0, 1.8 * gene_lane_count($gene_tsv)) : 0;
write_gnuplot(
    gp_file     => $gp_file,
    png_file    => $png_file,
    plot_tsv    => $plot_tsv,
    gene_tsv    => ($has_gtf ? $gene_tsv : ''),
    title       => ($opt{title} || "$opt{snp} local locus"),
    width       => $opt{width},
    height      => $opt{height},
    top_logp    => $opt{top_logp},
    labels      => \@labels,
    sig_y       => $sig_y,
    start       => $start,
    end         => $end,
    snp         => $opt{snp},
    target_chr  => $target_chr,
    target_bp   => $target_bp,
    target_markers => \@target_markers,
    ld_points    => $ld_points,
    ld_marker_symbol => lc($opt{ld_marker_symbol}),
    ld_marker_point_type => $ld_point_type,
    ld_marker_color => $opt{ld_marker_color},
    ld_display_mode => $opt{ld_display_mode},
    ld_population   => $opt{ld_population},
    ld_heatmap_colors => \@ld_heatmap_colors,
    ld_r2_points    => scalar(grep { exists $ld_r2_for{$_} && exists $found_ld_snp{$_} } keys %ld_r2_for),
    gene_height => $gene_height,
    use_zcolors => ($has_gtf && $has_zcols ? 1 : 0),
    colorbar_label => infer_effect_metric_label_from_cols(@resolved_zcols),
);

system($opt{gnuplot}, $gp_file) == 0
    or die "gnuplot failed for $gp_file\n";

open my $mf, '>', $manifest or die "Cannot write $manifest: $!\n";
print {$mf} join("\t", qw(METRIC VALUE)), "\n";
print {$mf} join("\t", 'cache_schema', 4), "\n";
print {$mf} join("\t", 'input', $opt{data}), "\n";
print {$mf} join("\t", 'png', $png_file), "\n";
print {$mf} join("\t", 'plot_tsv', $plot_tsv), "\n";
print {$mf} join("\t", 'gene_tsv', ($has_gtf ? $gene_tsv : '')), "\n";
print {$mf} join("\t", 'gtf', ($has_gtf ? $opt{gtf} : '')), "\n";
print {$mf} join("\t", 'has_gtf', ($has_gtf ? 1 : 0)), "\n";
print {$mf} join("\t", 'snp', $opt{snp}), "\n";
print {$mf} join("\t", 'label_snps', join(',', @label_snps)), "\n";
print {$mf} join("\t", 'label_snp_positions', join(',', map { $_->{snp} . ':' . $_->{bp} } @target_markers)), "\n";
print {$mf} join("\t", 'ld_snps', join(',', @ld_snps)), "\n";
print {$mf} join("\t", 'ld_snps_found', join(',', map { $found_ld_snp{$_} } sort keys %found_ld_snp)), "\n";
print {$mf} join("\t", 'ld_points_plotted', $ld_points), "\n";
print {$mf} join("\t", 'ld_marker_symbol', lc($opt{ld_marker_symbol})), "\n";
print {$mf} join("\t", 'ld_marker_color', $opt{ld_marker_color}), "\n";
print {$mf} join("\t", 'ld_display_mode', $opt{ld_display_mode}), "\n";
print {$mf} join("\t", 'ld_r2_values', ($opt{ld_r2_values} // '')), "\n";
print {$mf} join("\t", 'ld_population', $opt{ld_population}), "\n";
print {$mf} join("\t", 'ld_heatmap_colors', join(',', @ld_heatmap_colors)), "\n";
print {$mf} join("\t", 'chr', $target_chr), "\n";
print {$mf} join("\t", 'bp', $target_bp), "\n";
print {$mf} join("\t", 'window_bp', $window_bp), "\n";
print {$mf} join("\t", 'pcols', join(',', @pcols)), "\n";
print {$mf} join("\t", 'zcols', join(',', @zcols)), "\n";
print {$mf} join("\t", 'labels', ($opt{labels} // '')), "\n";
print {$mf} join("\t", 'title', ($opt{title} // '')), "\n";
print {$mf} join("\t", 'width', $opt{width}), "\n";
print {$mf} join("\t", 'height', $opt{height}), "\n";
print {$mf} join("\t", 'top_logp', $opt{top_logp}), "\n";
print {$mf} join("\t", 'sig', $opt{sig}), "\n";
print {$mf} join("\t", 'hide_y_axis', ($opt{hide_y_axis} ? 1 : 0)), "\n";
print {$mf} join("\t", 'bottom_snp_label', ($opt{bottom_snp_label} // '')), "\n";
print {$mf} join("\t", 'bottom_gene_label', ($opt{bottom_gene_label} // '')), "\n";
print {$mf} join("\t", 'rows_in_window', scalar(@locus)), "\n";
print {$mf} join("\t", 'points_plotted', $kept_points), "\n";
print {$mf} join("\t", 'gene_rows', $gene_rows), "\n";
close $mf or die "Cannot close $manifest: $!\n";

print "PNG\t$png_file\n";
print "MANIFEST\t$manifest\n";

sub write_gnuplot {
    my (%args) = @_;
    open my $gp, '>', $args{gp_file} or die "Cannot write $args{gp_file}: $!\n";
    my $track_count = scalar @{ $args{labels} };
    my $ymax = $track_count * $args{top_logp} + 1.2;
    my $ymin = $args{gene_tsv} ? -$args{gene_height} : (($args{bottom_gene_label} || $args{bottom_snp_label}) ? -6.8 : -2.6);
    my $chr_color = chromosome_color($args{target_chr});
    my $xmid = ($args{start} + $args{end}) / 2;
    print {$gp} "set terminal png enhanced size $args{width},$args{height}\n";
    print {$gp} "set output '" . escape_gp($args{png_file}) . "'\n";
    print {$gp} "set datafile separator '\\t'\n";
    print {$gp} "set title \"" . escape_gp($args{title}) . "\"\n";
    print {$gp} "set xlabel 'Position (bp)'\n";
    print {$gp} "set xrange [$args{start}:$args{end}]\n";
    print {$gp} "set yrange [$ymin:$ymax]\n";
    print {$gp} "set border " . ($args{hide_y_axis} ? 1 : 3) . "\n";
    print {$gp} "set lmargin " . ($args{hide_y_axis} ? 2 : 8) . "\n";
    print {$gp} "set rmargin 6\n";
    print {$gp} "set tics out nomirror\n";
    print {$gp} "set grid ytics lc rgb '#dddddd' dt 2\n";
    print {$gp} "unset key\n";
    my @ytics = repeated_panel_ytics_local($track_count, $args{top_logp});
    if ($args{hide_y_axis}) {
        print {$gp} "unset ylabel\n";
        print {$gp} "unset ytics\n";
        print {$gp} "unset mytics\n";
    }
    else {
        print {$gp} "set ylabel '-log10(P)'\n";
        print {$gp} "set ytics (" . join(', ', @ytics) . ")\n";
        print {$gp} "set mytics 2\n";
    }
    my $refline_ymax = ($args{gene_tsv} && -s $args{gene_tsv}) ? $ymax : 0;
    my $target_arrow_id = 1;
    for my $marker (@{ $args{target_markers} || [{ snp => $args{snp}, bp => $args{target_bp} }] }) {
        print {$gp} "set arrow $target_arrow_id from $marker->{bp},$ymin to $marker->{bp},$refline_ymax nohead lc rgb '#a0a0a0' dt 2 lw 1\n";
        $target_arrow_id++;
    }

    my $arrow_id = 10;
    for my $i (0 .. $#{ $args{labels} }) {
        my $base = $i * $args{top_logp};
        print {$gp} "set arrow $arrow_id from $args{start}," . ($base + $args{sig_y}) . " to $args{end}," . ($base + $args{sig_y}) . " nohead dt 2 lc rgb '#777777' lw 1\n";
        $arrow_id++;
        next unless $i > 0;
        print {$gp} "set arrow $arrow_id from $args{start},$base to $args{end},$base nohead lc rgb '#bbbbbb' lw 1\n";
        $arrow_id++;
    }
    my $draw_ld_markers = $args{ld_points}
        && $args{ld_display_mode} =~ /^(?:markers|both)$/;
    my $draw_ld_heatmap = $args{ld_r2_points}
        && $args{ld_display_mode} =~ /^(?:heatmap|both)$/;
    if ($draw_ld_markers) {
        print {$gp} "set label 900 '" . escape_gp("LD-linked SNP (" . $args{ld_marker_symbol} . ")") . "' at graph 0.015,0.985 left front tc rgb '" . escape_gp($args{ld_marker_color}) . "' font ',9'\n";
    }
    write_ld_heatmap_inset($gp, \%args) if $draw_ld_heatmap;
    my @ld_layers;
    if ($draw_ld_heatmap) {
        push @ld_layers, "'" . escape_gp($args{plot_tsv})
            . "' using ((\$9>=0)?\$1:1/0):2:10 with points pt 7 ps 1.35 lw 1 lc rgb variable";
    }
    if ($draw_ld_markers) {
        push @ld_layers, "'" . escape_gp($args{plot_tsv})
            . "' using ((\$8==1)?\$1:1/0):2 with points pt $args{ld_marker_point_type} ps 1.8 lw 2 lc rgb '"
            . escape_gp($args{ld_marker_color}) . "'";
    }
    for my $i (0 .. $#{ $args{labels} }) {
        my $base = $i * $args{top_logp};
        my $panel_y = $base + $args{top_logp} - 0.95;
        print {$gp} "set label " . (100 + $i) . " \"" . escape_gp($args{labels}[$i]) . "\" at $xmid,$panel_y center font ',18'\n";
    }

    if ($args{gene_tsv} && -s $args{gene_tsv}) {
        open my $gt, '<', $args{gene_tsv} or die "Cannot read $args{gene_tsv}: $!\n";
        <$gt>;
        my $obj = 1000;
        my $lab = 2000;
        my $min_exon_bp = max_num(600, int(($args{end} - $args{start}) * 0.003));
        my @gene_colors = (
            '#f4a6a6', '#f7c97f', '#d8d35f', '#9fd27a', '#71c9b8',
            '#86b6f6', '#a995e8', '#de9de6', '#f4a8c4', '#c8b08b',
            '#ef8b62', '#9ec3a5', '#7fb3d5', '#c39bd3', '#f8c471',
        );
        my %gene_color_for;
        my %lane_label_count;
        my @lane_label_offsets = (0.12, 0.30, 0.48, 0.66);
        while (my $line = <$gt>) {
            chomp $line;
            next unless length $line;
            my ($type, $gene, $gs, $ge, $lane) = split /\t/, $line, -1;
            my $track_top = -0.55 - (1.35 * $lane);
            my $line_top = $track_top - 0.12;
            my $line_bot = $line_top - 0.22;
            my $exon_top = $track_top + 0.02;
            my $exon_bot = $exon_top - 0.52;
            my $mid = ($gs + $ge) / 2;
            my $gene_color = $gene_color_for{$gene};
            if (!defined $gene_color) {
                my $color_idx = scalar(keys %gene_color_for) % @gene_colors;
                $gene_color = $gene_colors[$color_idx];
                $gene_color_for{$gene} = $gene_color;
            }
            if ($type eq 'gene') {
                my $lane_rank = $lane_label_count{$lane} // 0;
                my @below_label_offsets = (0.08, 0.22, 0.36, 0.50);
                my $label_y = $line_bot - $below_label_offsets[$lane_rank % @below_label_offsets];
                $lane_label_count{$lane} = $lane_rank + 1;
                print {$gp} "set object $obj rect from $gs,$line_bot to $ge,$line_top fc rgb '$gene_color' fillstyle solid 0.72 border lc rgb '$gene_color' lw 1 behind\n";
                print {$gp} "set label $lab \"" . escape_gp($gene) . "\" at $mid,$label_y center front tc rgb '#111111' font '" . italic_font_spec_gp(10) . "'\n";
                $lab++;
            }
            else {
                my $disp_start = $gs;
                my $disp_end = $ge;
                if (($disp_end - $disp_start) < $min_exon_bp) {
                    my $half_w = $min_exon_bp / 2;
                    $disp_start = $mid - $half_w;
                    $disp_end = $mid + $half_w;
                }
                print {$gp} "set object $obj rect from $disp_start,$exon_bot to $disp_end,$exon_top fc rgb '$gene_color' fillstyle solid 1.0 border lc rgb '$gene_color' lw 1 front\n";
            }
            $obj++;
        }
        close $gt;
    }

    if ($args{gene_tsv} && -s $args{gene_tsv}) {
        my @markers = sort { $a->{bp} <=> $b->{bp} } @{ $args{target_markers} || [{ snp => $args{snp}, bp => $args{target_bp} }] };
        my $plot_span = $args{end} - $args{start};
        my $edge_pad = $plot_span * 0.06;
        my $min_gap = $plot_span * 0.12;
        if (@markers > 1) {
            my $fit_gap = ($plot_span - 2 * $edge_pad) / (@markers - 1);
            $min_gap = $fit_gap if $fit_gap < $min_gap;
        }
        my @label_x = map { $_->{bp} } @markers;
        for my $i (1 .. $#label_x) {
            my $needed = $label_x[$i - 1] + $min_gap;
            $label_x[$i] = $needed if $label_x[$i] < $needed;
        }
        my $right_limit = $args{end} - $edge_pad;
        if (@label_x && $label_x[-1] > $right_limit) {
            $label_x[-1] = $right_limit;
            for (my $i = $#label_x - 1; $i >= 0; $i--) {
                my $needed = $label_x[$i + 1] - $min_gap;
                $label_x[$i] = $needed if $label_x[$i] > $needed;
            }
        }
        my $left_limit = $args{start} + $edge_pad;
        if (@label_x && $label_x[0] < $left_limit) {
            my $shift = $left_limit - $label_x[0];
            $_ += $shift for @label_x;
        }
        my $label_y = $ymax - 0.22;
        my $leader_y = $ymax - 0.85;
        for my $i (0 .. $#markers) {
            my $label_id = 1 + $i;
            my $leader_id = 500 + $i;
            print {$gp} "set arrow $leader_id from $markers[$i]{bp},$leader_y to $label_x[$i]," . ($label_y - 0.12) . " nohead lc rgb '#777777' dt 3 lw 1 front\n";
            print {$gp} "set label $label_id \"" . escape_gp($markers[$i]{snp}) . "\" at $label_x[$i],$label_y center font ',10' front\n";
        }
        print {$gp} "set xlabel 'Chromosome " . escape_gp($args{target_chr}) . "'\n";
        if ($args{use_zcolors}) {
            print {$gp} "set cbrange [-8:8]\n";
            print {$gp} "set cbtics ('-8' -8, '0' 0, '8' 8)\n";
            print {$gp} "set cblabel '" . escape_gp($args{colorbar_label} || 'Effect metric') . "'\n";
            print {$gp} "set colorbox vertical user origin 0.94,0.12 size 0.02,0.76\n";
            print {$gp} "set palette defined (-8 '#63d67f', -4 '#63d8d2', 0 '#ffbf00', 4 '#ff5b00', 8 '#df1f2d')\n";
            my @plots = ("'" . escape_gp($args{plot_tsv}) . "' using 1:2:7 with points pt 7 ps 0.72 lc palette");
            push @plots, @ld_layers;
            print {$gp} "plot " . join(', ', @plots) . "\n";
        }
        else {
            my @palette = sas_chr_palette();
            print {$gp} "set palette maxcolors " . scalar(@palette) . " defined (" .
                join(', ', map { sprintf('%d "%s"', $_ + 1, $palette[$_]) } 0 .. $#palette) . ")\n";
            print {$gp} "unset colorbox\n";
            my @plots = ("'" . escape_gp($args{plot_tsv}) . "' using 1:2:(\$3+1) with points pt 7 ps 0.72 lc palette");
            push @plots, @ld_layers;
            print {$gp} "plot " . join(', ', @plots) . "\n";
        }
    }
    else {
        print {$gp} "unset xlabel\n";
        print {$gp} "unset xtics\n";
        my $snp_label = $args{bottom_snp_label} || $args{snp};
        my $gene_label = $args{bottom_gene_label} || '';
        my $label_dx = ($args{end} - $args{start}) * 0.016;
        my $label_y = $ymin + (length($gene_label) ? 4.75 : 2.35);
        print {$gp} "set label 1 \"" . escape_gp($snp_label) . "\" at " . ($args{target_bp} - $label_dx) . ",$label_y center rotate by 90 font ',11'\n";
        if (length $gene_label) {
            print {$gp} "set label 2 \"" . escape_gp($gene_label) . "\" at " . ($args{target_bp} + $label_dx) . ",$label_y center rotate by 90 font '" . italic_font_spec_gp(11) . "'\n";
        }
        print {$gp} "unset colorbox\n";
        my @plots = ("'" . escape_gp($args{plot_tsv}) . "' using 1:2 with points pt 7 ps 0.9 lc rgb '" . escape_gp($chr_color) . "'");
        push @plots, @ld_layers;
        print {$gp} "plot " . join(', ', @plots) . "\n";
    }
    close $gp or die "Cannot close $args{gp_file}: $!\n";
}

sub write_ld_heatmap_inset {
    my ($gp, $args) = @_;
    my @colors = @{ $args->{ld_heatmap_colors} || [] };
    return unless @colors >= 2;
    my $segments = 8;
    my ($x0, $x1, $y0, $y1) = (0.70, 0.90, 0.935, 0.957);
    my $dx = ($x1 - $x0) / $segments;
    for my $i (0 .. $segments - 1) {
        my $r2 = ($i + 0.5) / $segments;
        my $color = ld_color_for_r2($r2, \@colors);
        my $left = $x0 + $i * $dx;
        my $right = $left + $dx;
        my $id = 8000 + $i;
        print {$gp} "set object $id rect from graph $left,$y0 to graph $right,$y1 fc rgb '$color' fillstyle solid 1.0 border lc rgb '$color' front\n";
    }
    my $title = 'LD r^2 (' . ($args->{ld_population} || 'EUR') . ')';
    print {$gp} "set label 8900 '" . escape_gp($title) . "' at graph " . (($x0 + $x1) / 2) . "," . ($y1 + 0.022) . " center front font ',9' tc rgb '#222222'\n";
    print {$gp} "set label 8901 '0' at graph $x0," . ($y0 - 0.010) . " center front font ',8' tc rgb '#222222'\n";
    print {$gp} "set label 8902 '0.5' at graph " . (($x0 + $x1) / 2) . "," . ($y0 - 0.010) . " center front font ',8' tc rgb '#222222'\n";
    print {$gp} "set label 8903 '1' at graph $x1," . ($y0 - 0.010) . " center front font ',8' tc rgb '#222222'\n";
}

sub allocate_lane {
    my ($lanes, $start, $end) = @_;
    for my $i (0 .. $#{$lanes}) {
        next if $start <= $lanes->[$i];
        $lanes->[$i] = $end;
        return $i;
    }
    push @{$lanes}, $end;
    return $#{$lanes};
}

sub lane_reservation_range {
    my (%args) = @_;
    my $gene = $args{gene} // 'GENE';
    my $gene_start = $args{gene_start};
    my $gene_end = $args{gene_end};
    my $plot_start = $args{plot_start};
    my $plot_end = $args{plot_end};
    my $plot_width = $args{plot_width} || 1500;
    my $plot_bp_span = $args{plot_bp_span} || max_num(1, $plot_end - $plot_start + 1);

    my $gene_span = max_num(1, $gene_end - $gene_start);
    my $chars = length($gene);
    my $usable_width_px = max_num(600, $plot_width - 180);
    my $px_per_char = 8.5;
    my $label_px = max_num(40, int($chars * $px_per_char));
    my $label_bp_span = int(($label_px / $usable_width_px) * $plot_bp_span);
    my $clearance_bp = max_num(
        int($plot_bp_span * 0.010),
        max_num(int($gene_span * 0.35), 1200),
    );
    my $reserve_span = max_num($gene_span + 2 * $clearance_bp, $label_bp_span + 2 * $clearance_bp);
    my $mid = ($gene_start + $gene_end) / 2;
    my $lane_start = int($mid - $reserve_span / 2);
    my $lane_end = int($mid + $reserve_span / 2);
    $lane_start = $plot_start if $lane_start < $plot_start;
    $lane_end = $plot_end if $lane_end > $plot_end;
    if ($lane_end <= $lane_start) {
        $lane_end = $lane_start + 1;
    }
    return ($lane_start, $lane_end);
}

sub gene_lane_count {
    my ($path) = @_;
    return 0 unless $path && -s $path;
    open my $fh, '<', $path or return 0;
    <$fh>;
    my $max_lane = -1;
    while (my $line = <$fh>) {
        chomp $line;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        my $lane = numeric($f[4]);
        next unless defined $lane;
        $max_lane = $lane if $lane > $max_lane;
    }
    close $fh;
    return $max_lane + 1;
}

sub repeated_panel_ytics_local {
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

sub max_num {
    my ($a, $b) = @_;
    return $a >= $b ? $a : $b;
}

sub cap_num {
    my ($x, $minv, $maxv) = @_;
    return $minv if $x < $minv;
    return $maxv if $x > $maxv;
    return $x;
}

sub ld_rgb_integer {
    my ($r2, $colors) = @_;
    my $hex = ld_color_for_r2($r2, $colors);
    $hex =~ s/^#//;
    return hex($hex);
}

sub ld_color_for_r2 {
    my ($r2, $colors) = @_;
    $r2 = cap_num($r2, 0, 1);
    my $scaled = $r2 * (@{$colors} - 1);
    my $left = int($scaled);
    $left = @{$colors} - 2 if $left >= @{$colors} - 1;
    my $fraction = $scaled - $left;
    $fraction = 1 if $r2 >= 1;
    my @a = hex_color_rgb($colors->[$left]);
    my @b = hex_color_rgb($colors->[$left + 1]);
    my @rgb = map { int($a[$_] + ($b[$_] - $a[$_]) * $fraction + 0.5) } 0 .. 2;
    return sprintf('#%02x%02x%02x', @rgb);
}

sub hex_color_rgb {
    my ($hex) = @_;
    $hex =~ s/^#//;
    return map { hex($_) } ($hex =~ /(..)(..)(..)/);
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

sub infer_effect_metric_label_from_cols {
    my (@cols) = @_;
    @cols = grep { defined $_ && length $_ } @cols;
    return 'Effect metric' unless @cols;
    my $all_z = 1;
    my $all_beta = 1;
    my $all_or = 1;
    for my $col (@cols) {
        my $u = uc($col);
        $all_z = 0 unless $u =~ /(?:^|_)Z(?:SCORE)?(?:_|$)/;
        $all_beta = 0 unless $u =~ /BETA/;
        $all_or = 0 unless $u =~ /(?:^|_)OR(?:_|$)|ODDSRATIO/;
    }
    return 'Z score' if $all_z;
    return 'Beta' if $all_beta;
    return 'Odds ratio' if $all_or;
    return 'Effect metric';
}

sub chromosome_color {
    my ($chr) = @_;
    my @palette = genomewide_palette();
    my $ord = chr_palette_index($chr);
    return $palette[($ord - 1) % @palette];
}

sub chr_palette_index {
    my ($chr) = @_;
    $chr = normalize_chr($chr);
    return 23 if $chr eq 'X';
    return $chr if length($chr) && $chr =~ /^\d+$/ && $chr > 0;
    return 1;
}

sub genomewide_palette {
    return (
        sas_chr_palette()
    );
}

sub sas_chr_palette {
    return (
        '#0072bd', '#d95319', '#edb120', '#7e2f8e',
        '#77ac30', '#4dbeee', '#a2142f',
    );
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

sub extract_requested_numeric {
    my ($requested, $row, $idx) = @_;
    if (exists $idx->{$requested}) {
        return numeric($row->[ $idx->{$requested} ]);
    }
    my @cand = alias_candidates($requested);
    for my $cand (@cand) {
        return numeric($row->[ $idx->{$cand} ]) if exists $idx->{$cand};
    }
    if ($requested =~ /^(.*)_FEMALE_Z$/) {
        return beta_se_to_z($1 . '_GROUP1_BETA', $1 . '_GROUP1_SE', $row, $idx);
    }
    if ($requested =~ /^(.*)_MALE_Z$/) {
        return beta_se_to_z($1 . '_GROUP2_BETA', $1 . '_GROUP2_SE', $row, $idx);
    }
    if ($requested =~ /^(.*)_STD_Z$/) {
        return beta_se_to_z($1 . '_DIFF_BETA', $1 . '_DIFF_SE', $row, $idx);
    }
    return undef;
}

sub beta_se_to_z {
    my ($beta_col, $se_col, $row, $idx) = @_;
    return undef unless exists $idx->{$beta_col} && exists $idx->{$se_col};
    my $beta = numeric($row->[ $idx->{$beta_col} ]);
    my $se = numeric($row->[ $idx->{$se_col} ]);
    return undef unless defined $beta && defined $se && $se > 0;
    return $beta / $se;
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

sub safe_neglog10 {
    my ($value) = @_;
    my $num = numeric($value);
    return undef unless defined $num && $num > 0;
    return -log($num) / log(10);
}

sub ld_marker_point_type {
    my ($name) = @_;
    $name = lc(trim($name || 'star'));
    my %point_type = (
        plus     => 1,
        cross    => 2,
        star     => 3,
        square   => 5,
        circle   => 7,
        triangle => 9,
        diamond  => 13,
    );
    die "Unsupported --ld-marker-symbol '$name'; use star, plus, cross, circle, square, triangle, or diamond\n"
        unless exists $point_type{$name};
    return $point_type{$name};
}

sub ld_marker_legend_text {
    my ($name) = @_;
    $name = lc($name || 'star');
    return '+' if $name eq 'plus';
    return 'x' if $name eq 'cross';
    return 'o' if $name eq 'circle';
    return 's' if $name eq 'square';
    return '^' if $name eq 'triangle';
    return 'd' if $name eq 'diamond';
    return '*';
}

sub nice_label {
    my ($text) = @_;
    $text =~ s/_/ /g;
    $text =~ s/\bSTD\b/standardized/ig;
    $text =~ s/\bDIFF\b/diff/ig;
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

sub italic_font_spec_gp {
    my ($size) = @_;
    $size ||= 10;
    my @candidates = (
        'C:/Windows/Fonts/timesi.ttf',
        'C:/Windows/Fonts/ariali.ttf',
    );
    for my $cand (@candidates) {
        return $cand . "," . $size if -f $cand;
    }
    return "," . $size;
}
