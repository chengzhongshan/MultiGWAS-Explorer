#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use File::Basename qw(basename);
use File::Path qw(make_path);
use Text::ParseWords qw(parse_line);

sub usage {
    return <<"USAGE";
Usage:
  perl pdl_gunplot_forest.pl --csv top_hits.csv --out-prefix prefix [options]

Options:
  --track-ids A|B|C         Required.
  --track-labels A|B|C      Required.
  --track-beta-vars A|B|C   Required.
  --track-se-vars A|B|C     Required.
  --track-p-vars A|B|C      Required.
  --title TEXT              Optional title.
  --width N                 Default: 1100
  --height N                Optional fixed panel height.
  --dotsize N               Default: 8
  --y-font-size N           Default: 12
  --min-axis FLOAT          Default: 0.4
  --max-axis FLOAT          Default: 1.6
  --xaxis-ticks TEXT        Example: 0.4 to 1.6 by 0.2
  --default-hit-class TEXT  Default: DIFFERENTIAL
  --gnuplot PATH            Default: gnuplot
USAGE
}

my %opt = (
    width => 1100,
    height => '',
    dotsize => 8,
    y_font_size => 12,
    min_axis => 0.4,
    max_axis => 1.6,
    xaxis_ticks => '0.4 to 1.6 by 0.2',
    default_hit_class => 'DIFFERENTIAL',
    gnuplot => 'gnuplot',
    title => 'Forest plot',
);

GetOptions(
    'csv=s'               => \$opt{csv},
    'out-prefix=s'        => \$opt{out_prefix},
    'track-ids=s'         => \$opt{track_ids},
    'track-labels=s'      => \$opt{track_labels},
    'track-beta-vars=s'   => \$opt{track_beta_vars},
    'track-se-vars=s'     => \$opt{track_se_vars},
    'track-p-vars=s'      => \$opt{track_p_vars},
    'title=s'             => \$opt{title},
    'width=i'             => \$opt{width},
    'height=s'            => \$opt{height},
    'dotsize=f'           => \$opt{dotsize},
    'y-font-size=i'       => \$opt{y_font_size},
    'min-axis=f'          => \$opt{min_axis},
    'max-axis=f'          => \$opt{max_axis},
    'xaxis-ticks=s'       => \$opt{xaxis_ticks},
    'default-hit-class=s' => \$opt{default_hit_class},
    'gnuplot=s'           => \$opt{gnuplot},
) or die usage();

die usage() unless $opt{csv} && $opt{out_prefix} && $opt{track_ids} && $opt{track_labels}
    && $opt{track_beta_vars} && $opt{track_se_vars} && $opt{track_p_vars};
die "Input CSV not found: $opt{csv}\n" unless -s $opt{csv};

my @track_ids = split /\|/, $opt{track_ids};
my @track_labels = split /\|/, $opt{track_labels};
my @track_beta_vars = split /\|/, $opt{track_beta_vars};
my @track_se_vars = split /\|/, $opt{track_se_vars};
my @track_p_vars = split /\|/, $opt{track_p_vars};
my $track_count = scalar(@track_ids);
die "Forest track arrays are inconsistent\n"
    unless $track_count
        && @track_labels == $track_count
        && @track_beta_vars == $track_count
        && @track_se_vars == $track_count
        && @track_p_vars == $track_count;

my ($header, $rows) = read_csv_rows($opt{csv});
my %header_index = map { $header->[$_] => $_ } 0 .. $#$header;
die "Forest CSV is missing SNP column\n" unless exists $header_index{SNP};

my %snp_seen;
my @row_order;
for my $row (@{$rows}) {
    my $snp = trim($row->{SNP} // '');
    next unless length $snp;
    next if $snp_seen{$snp}++;
    push @row_order, $row;
}
die "No SNP rows were loaded from $opt{csv}\n" unless @row_order;

my $single_snp_mode = @row_order == 1 ? 1 : 0;
my $manifest_path = $opt{out_prefix} . '.manifest.tsv';

open my $mf, '>', $manifest_path or die "Cannot write $manifest_path: $!\n";
print {$mf} join("\t", qw(track_order track_id track_label png_file)), "\n";

if ($single_snp_mode) {
    my $row = $row_order[0];
    my @points;
    for my $i (0 .. $#track_ids) {
        my $beta = numeric($row->{ $track_beta_vars[$i] });
        my $se   = numeric($row->{ $track_se_vars[$i] });
        my $p    = numeric($row->{ $track_p_vars[$i] });
        next unless defined $beta && defined $se && $se > 0 && defined $p && $p > 0;
        push @points, build_point(
            left_label    => $track_labels[$i],
            right_label   => '',
            beta          => $beta,
            se            => $se,
            p             => $p,
            hit_class     => $row->{hit_class} || $opt{default_hit_class},
        );
    }
    die "No non-missing cohort rows were available for the single-SNP forest plot\n" unless @points;
    my $single_title = single_snp_title($row, $opt{title});
    my $png_path = $opt{out_prefix} . '_single_snp.png';
    render_panel(
        panel_id     => 'single_snp',
        panel_title  => $single_title,
        png_path     => $png_path,
        points       => \@points,
        width        => $opt{width},
        height       => effective_panel_height($opt{height}, scalar(@points), 1),
        dotsize      => $opt{dotsize},
        y_font_size  => $opt{y_font_size},
        min_axis     => $opt{min_axis},
        max_axis     => $opt{max_axis},
        xaxis_ticks  => $opt{xaxis_ticks},
        gnuplot      => $opt{gnuplot},
        single_snp_mode => 1,
    );
    print {$mf} join("\t", 1, 'single_snp', $single_title, basename($png_path)), "\n";
} else {
    my @prepared_rows = sort {
        class_sort_key($a->{hit_class} || $opt{default_hit_class}) <=> class_sort_key($b->{hit_class} || $opt{default_hit_class})
            ||
        numeric_or_large($a->{hit_order}) <=> numeric_or_large($b->{hit_order})
            ||
        (($a->{SNP} || '') cmp ($b->{SNP} || ''))
    } @row_order;

    for my $track_i (0 .. $#track_ids) {
        my @points;
        my @separator_values;
        my $prev_class = '';
        for my $row (@prepared_rows) {
            my $beta = numeric($row->{ $track_beta_vars[$track_i] });
            my $se   = numeric($row->{ $track_se_vars[$track_i] });
            my $p    = numeric($row->{ $track_p_vars[$track_i] });
            next unless defined $beta && defined $se && $se > 0 && defined $p && $p > 0;
            my $class = $row->{hit_class} || $opt{default_hit_class};
            if (@points && uc(trim($class)) ne uc(trim($prev_class))) {
                push @separator_values, scalar(@points) + 0.5;
            }
            push @points, build_point(
                left_label    => trim($row->{SNP} // ''),
                right_label   => normalize_gene_label($row->{gene}),
                beta          => $beta,
                se            => $se,
                p             => $p,
                hit_class     => $class,
            );
            $prev_class = $class;
        }
        next unless @points;
        my $track_id = safe_name($track_ids[$track_i]);
        my $track_label = $track_labels[$track_i];
        my $png_path = $opt{out_prefix} . '_' . $track_id . '.png';
        render_panel(
            panel_id     => $track_id,
            panel_title  => $track_label,
            png_path     => $png_path,
            points       => \@points,
            width        => $opt{width},
            height       => effective_panel_height($opt{height}, scalar(@points), 0),
            dotsize      => $opt{dotsize},
            y_font_size  => $opt{y_font_size},
            min_axis     => $opt{min_axis},
            max_axis     => $opt{max_axis},
            xaxis_ticks  => $opt{xaxis_ticks},
            gnuplot      => $opt{gnuplot},
            separator_values => \@separator_values,
            single_snp_mode => 0,
        );
        print {$mf} join("\t", $track_i + 1, $track_id, $track_label, basename($png_path)), "\n";
    }
}
close $mf or die "Cannot close $manifest_path: $!\n";

print "Wrote forest manifest: $manifest_path\n";

sub render_panel {
    my (%args) = @_;
    my $png_path = $args{png_path};
    my $data_path = $png_path;
    $data_path =~ s/\.png$/.data.tsv/i;
    my $gp_path = $png_path;
    $gp_path =~ s/\.png$/.gp/i;
    my $gp_png_path = gnuplot_path($png_path);
    my $gp_data_path = gnuplot_path($data_path);

    open my $df, '>', $data_path or die "Cannot write $data_path: $!\n";
    print {$df} join("\t", qw(Y OR LOW HIGH SIG CLASS STAR LEFT_LABEL RIGHT_LABEL)), "\n";
    for my $i (0 .. $#{ $args{points} || [] }) {
        my $point = $args{points}[$i];
        my $y = $i + 1;
        print {$df} join(
            "\t",
            $y,
            map { defined $_ ? $_ : '' }
            (
                sprintf('%.8f', $point->{or}),
                sprintf('%.8f', $point->{low}),
                sprintf('%.8f', $point->{high}),
                $point->{sig},
                $point->{class},
                $point->{star},
                $point->{left_label},
                $point->{right_label},
            )
        ), "\n";
    }
    close $df or die "Cannot close $data_path: $!\n";

    my $n = scalar(@{ $args{points} || [] });
    my @ytics;
    my @y2tics;
    for my $i (0 .. $#{ $args{points} || [] }) {
        my $point = $args{points}[$i];
        push @ytics, sprintf('"%s" %d', gp_escape($point->{left_label}), $i + 1);
        if (!$args{single_snp_mode}) {
            push @y2tics, sprintf('"%s" %d', gp_enhanced_italic($point->{right_label}), $i + 1);
        }
    }
    my $ytics_text = join(", ", @ytics);
    my $y2tics_text = join(", ", @y2tics);
    my $xtics_cmd = build_xtics_cmd($args{xaxis_ticks});
    my $point_size = point_size_from_dotsize($args{dotsize});
    my @plot_clauses;
    if ($args{single_snp_mode}) {
        @plot_clauses = (
            qq{"$data_path" using 2:1:3:4 with xerrorbars lw 1.2 lc rgb "#4f67b0"},
            qq{"" using 2:1 with points pt 7 ps $point_size lc rgb "#4f67b0"},
            qq{"" using ((\$5 == 1) ? \$2 : 1/0):1:7 with labels center tc rgb "#111111" font "Arial,@{[$args{y_font_size} + 1]}"},
        );
    } else {
        @plot_clauses = (
            qq{"$data_path" using ((strcol(6) eq "COMMON") ? \$2 : 1/0):1:3:4 with xerrorbars lw 1.2 lc rgb "#4f67b0"},
            qq{"" using ((strcol(6) eq "COMMON") ? \$2 : 1/0):1 with points pt 7 ps $point_size lc rgb "#4f67b0"},
            qq{"" using ((strcol(6) eq "DIFFERENTIAL") ? \$2 : 1/0):1:3:4 with xerrorbars lw 1.2 lc rgb "#c0504d"},
            qq{"" using ((strcol(6) eq "DIFFERENTIAL") ? \$2 : 1/0):1 with points pt 7 ps $point_size lc rgb "#c0504d"},
            qq{"" using ((strcol(6) ne "COMMON" && strcol(6) ne "DIFFERENTIAL") ? \$2 : 1/0):1:3:4 with xerrorbars lw 1.2 lc rgb "#666666"},
            qq{"" using ((strcol(6) ne "COMMON" && strcol(6) ne "DIFFERENTIAL") ? \$2 : 1/0):1 with points pt 7 ps $point_size lc rgb "#666666"},
            qq{"" using ((\$5 == 1) ? \$2 : 1/0):1:7 with labels center tc rgb "#111111" font "Arial,@{[$args{y_font_size} + 1]}"},
        );
    }
    my @separator_cmds;
    for my $idx (0 .. $#{ $args{separator_values} || [] }) {
        my $y = $args{separator_values}[$idx];
        push @separator_cmds,
            sprintf('set arrow %d from %s,%s to %s,%s nohead dt 2 lw 1 lc rgb "#9a9a9a"', $idx + 20, $args{min_axis}, $y, $args{max_axis}, $y);
    }

    open my $gp, '>', $gp_path or die "Cannot write $gp_path: $!\n";
    print {$gp} <<"GP";
set terminal pngcairo size $args{width},$args{height} enhanced font "Arial,$args{y_font_size}"
set output "$gp_png_path"
unset key
set title "${\gp_escape($args{panel_title})}"
set xrange [$args{min_axis}:$args{max_axis}]
set yrange [0.5:$n+0.5]
set ytics nomirror font "Arial,$args{y_font_size}" ($ytics_text)
set xlabel "OR and 95% CI" font "Arial,$args{y_font_size}"
set ylabel "" font "Arial,$args{y_font_size}"
set border 3
set tmargin 2
set lmargin 18
set rmargin @{[$args{single_snp_mode} ? 4 : 18]}
set bmargin 4
unset grid
set xtics nomirror font "Arial,$args{y_font_size}"
@{[$args{single_snp_mode} ? '' : 'set y2range [0.5:'.$n.'+0.5]']}
@{[$args{single_snp_mode} ? '' : 'set y2tics nomirror font "Arial,'.$args{y_font_size}.'" textcolor rgb "#444444" ('.$y2tics_text.')']}
@{[$args{single_snp_mode} ? '' : 'set y2label "" font "Arial,'.$args{y_font_size}.'"']}
@{[$args{single_snp_mode} ? '' : 'set grid ytics lc rgb "#d9d9d9" dt 3 lw 1']}
set arrow 1 from 1, graph 0 to 1, graph 1 nohead lw 1 lc rgb "#777777"
$xtics_cmd
GP
    if (@separator_cmds) {
        print {$gp} join("\n", @separator_cmds), "\n";
    }
    for my $clause (@plot_clauses) {
        $clause =~ s/\Q$data_path\E/$gp_data_path/g;
    }
    print {$gp} "plot \\\n  " . join(",\\\n  ", @plot_clauses) . "\n";
    close $gp or die "Cannot close $gp_path: $!\n";
    system($args{gnuplot}, $gp_path) == 0
        or die "gnuplot failed for forest panel $gp_path\n";
}

sub build_point {
    my (%args) = @_;
    my $or = exp($args{beta});
    my $low = exp($args{beta} - 1.96 * $args{se});
    my $high = exp($args{beta} + 1.96 * $args{se});
    my $class = uc(trim($args{hit_class} // ''));
    $class = 'OTHER' unless length $class;
    my $sig = (defined($args{p}) && $args{p} > 0 && $args{p} < 5e-8) ? 1 : 0;
    return {
        left_label  => $args{left_label},
        right_label => $args{right_label},
        or          => $or,
        low         => $low,
        high        => $high,
        sig         => $sig,
        class       => $class,
        star        => ($sig ? '*' : ''),
    };
}

sub single_snp_title {
    my ($row, $prefix) = @_;
    my $snp = trim($row->{SNP} // '');
    my $gene = trim($row->{gene} // '');
    my $title = $snp;
    $title .= " ($gene)" if length($gene) && uc($gene) ne 'NA';
    return $title;
}

sub multi_snp_label {
    my ($row) = @_;
    my $snp = trim($row->{SNP} // '');
    my $gene = trim($row->{gene} // '');
    return $snp unless length($gene) && uc($gene) ne 'NA';
    return "$snp ($gene)";
}

sub normalize_gene_label {
    my ($gene) = @_;
    $gene = trim($gene // '');
    return 'NA' unless length($gene);
    return $gene;
}

sub effective_panel_height {
    my ($fixed_height, $n_rows, $single_snp_mode) = @_;
    return int($fixed_height) if defined($fixed_height) && length($fixed_height) && $fixed_height =~ /^\d+(?:\.\d+)?$/;
    return 320 + ($n_rows * 55) if $single_snp_mode;
    my $height = 240 + ($n_rows * 34);
    $height = 700 if $height < 700;
    return int($height);
}

sub class_sort_key {
    my ($class) = @_;
    $class = uc(trim($class // ''));
    return 1 if $class eq 'DIFFERENTIAL';
    return 2 if $class eq 'COMMON';
    return 3;
}

sub read_csv_rows {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    my $header = <$fh>;
    die "CSV is empty: $path\n" unless defined $header;
    chomp $header;
    $header =~ s/\r$//;
    my @cols = map { trim($_) } parse_line(',', 0, $header);
    my @rows;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @vals = parse_line(',', 0, $line);
        my %row;
        @row{@cols} = @vals;
        push @rows, \%row;
    }
    close $fh;
    return (\@cols, \@rows);
}

sub build_xtics_cmd {
    my ($text) = @_;
    $text //= '';
    if ($text =~ /^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+to\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+by\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s*$/i) {
        return sprintf('set xtics %s,%s,%s', $1, $3, $2);
    }
    return 'set xtics';
}

sub numeric {
    my ($value) = @_;
    return undef unless defined $value;
    $value = trim($value);
    return undef unless $value =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
    return 0 + $value;
}

sub numeric_or_large {
    my ($value) = @_;
    my $num = numeric($value);
    return defined($num) ? $num : 9e99;
}

sub gp_escape {
    my ($text) = @_;
    $text //= '';
    $text =~ s/\\/\\\\/g;
    $text =~ s/"/\\"/g;
    return $text;
}

sub gp_enhanced_italic {
    my ($text) = @_;
    $text = normalize_gene_label($text);
    $text =~ s/\\/\\\\/g;
    $text =~ s/"/\\"/g;
    $text =~ s/_/\\_/g;
    $text =~ s/\{/\\{/g;
    $text =~ s/\}/\\}/g;
    return "{/:Italic $text}";
}

sub point_size_from_dotsize {
    my ($dotsize) = @_;
    $dotsize = 8 unless defined $dotsize && $dotsize =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)$/;
    my $size = $dotsize / 5.5;
    $size = 0.8 if $size < 0.8;
    return sprintf('%.2f', $size);
}

sub safe_name {
    my ($text) = @_;
    $text //= '';
    $text =~ s/[^A-Za-z0-9._-]+/_/g;
    $text =~ s/^_+|_+$//g;
    return length($text) ? $text : 'item';
}

sub gnuplot_path {
    my ($path) = @_;
    $path //= '';
    $path =~ s{\\}{/}g;
    return $path;
}

sub trim {
    my ($text) = @_;
    $text = '' unless defined $text;
    $text =~ s/^\s+//;
    $text =~ s/\s+$//;
    return $text;
}
