#!/usr/bin/env perl
package LocalGTFLabelPlanner;
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use JSON::PP qw(encode_json);

# SAS text widths depend on the installed font. These conservative character
# widths, plus a gap, decide whether horizontal labels are safe before upload.
sub text_width_px {
    my ($label, $font_pt) = @_;
    my $em = 0;
    for my $char (split //, $label) {
        $em += $char =~ /[WM@%]/ ? 0.95
             : $char =~ /[ilI1.,:;|]/ ? 0.40
             : $char =~ /[A-Z]/ ? 0.72
             : $char =~ /[0-9]/ ? 0.62 : 0.62;
    }
    return ($em * $font_pt * 96 / 72) * 1.12 + 4;
}

sub place_centers {
    my ($items, $half_widths, $left, $right, $gap) = @_;
    my $count = @$items;
    my $needed = 2 * _sum(@$half_widths) + ($count - 1) * $gap;
    return undef if $needed > $right - $left;
    # Weighted isotonic projection: after subtracting each required gap,
    # nonoverlapping centers become a nondecreasing sequence. Pool adjacent
    # violations to share displacement instead of pushing only right labels.
    my @separation = (0);
    for my $i (1 .. $count - 1) {
        $separation[$i] = $separation[$i - 1] + $half_widths->[$i - 1]
            + $half_widths->[$i] + $gap;
    }
    my @blocks;
    for my $i (0 .. $count - 1) {
        push @blocks, { first => $i, last => $i,
            sum => $items->[$i]{x} - $separation[$i], count => 1 };
        while (@blocks >= 2
            && $blocks[-2]{sum} / $blocks[-2]{count}
                > $blocks[-1]{sum} / $blocks[-1]{count}) {
            my $right_block = pop @blocks;
            $blocks[-1]{last} = $right_block->{last};
            $blocks[-1]{sum} += $right_block->{sum};
            $blocks[-1]{count} += $right_block->{count};
        }
    }
    my $minimum = $left + $half_widths->[0];
    my $maximum = $right - $half_widths->[-1] - $separation[-1];
    my @center;
    for my $block (@blocks) {
        my $base = $block->{sum} / $block->{count};
        $base = $minimum if $base < $minimum;
        $base = $maximum if $base > $maximum;
        for my $i ($block->{first} .. $block->{last}) {
            $center[$i] = $base + $separation[$i];
        }
    }
    return \@center;
}

sub _sum { my $total = 0; $total += $_ for @_; return $total }

sub plan {
    my (%arg) = @_;
    my $targets = $arg{targets} || [];
    die "targets must be an array\n" unless ref($targets) eq 'ARRAY' && @$targets;
    my $layout = lc($arg{layout} || 'auto');
    die "layout must be auto, horizontal, or vertical\n"
        unless $layout =~ /^(?:auto|horizontal|vertical)$/;
    my $width = 0 + ($arg{design_width} || 950);
    my $height = 0 + ($arg{design_height} || 1000);
    my $font = 0 + ($arg{font_size} || 10);
    my $reference_bp = 0 + ($arg{reference_bp} || 0);
    my $window = 0 + ($arg{window_bp} || 0);
    die "design_width must be at least 300\n" if $width < 300;
    die "design_height must be at least 300\n" if $height < 300;
    die "font_size must be between 5 and 24\n" if $font < 5 || $font > 24;
    die "window_bp must be positive\n" if $window <= 0;
    my $x_min = $reference_bp - $window;
    $x_min = 1 if $x_min < 1;
    my $x_max = $reference_bp + $window;
    if ($arg{chromosome_length} && $arg{chromosome_length} < $x_max) {
        $x_max = 0 + $arg{chromosome_length};
    }
    die "invalid locus bounds\n" if $x_max <= $x_min;
    my $span = $x_max - $x_min;
    my $left = 0.11 * $width;
    my $right = 0.89 * $width;
    my $plot_width = $right - $left;
    my @items;
    for my $index (0 .. $#$targets) {
        my $target = $targets->[$index];
        my $label = $target->{snp} // '';
        die "invalid SNP label: $label\n" unless $label =~ /^[A-Za-z0-9_.:-]+$/;
        my $bp = $target->{bp};
        unless (defined($bp) && $bp =~ /^\d+(?:\.\d+)?$/
            && $bp >= $x_min && $bp <= $x_max) {
            return { layout => ($layout eq 'vertical' || @$targets > 1
                    && $layout ne 'horizontal' ? 'vertical' : 'horizontal'),
                positions => '', font_size => $font,
                headroom_frac => '', center_offset => '',
                reason => 'missing_or_outside_coordinate' };
        }
        push @items, { snp => $label, bp => 0 + $bp, order => $index,
            x => $left + ($bp - $x_min) / $span * $plot_width };
    }
    @items = sort { $a->{bp} <=> $b->{bp} || $a->{order} <=> $b->{order} } @items;
    my @horizontal_half = map { text_width_px($_->{snp}, $font) / 2 } @items;
    my $horizontal_gap = $font * 96 / 72 * 0.7 + 6;
    my $horizontal = place_centers(\@items, \@horizontal_half, $left, $right,
        $horizontal_gap);
    my $max_horizontal_shift = 0;
    if ($horizontal) {
        for my $i (0 .. $#items) {
            my $shift = abs($horizontal->[$i] - $items[$i]{x});
            $max_horizontal_shift = $shift if $shift > $max_horizontal_shift;
        }
    }
    # As in a minimum-separation space adjustment, keep moving overlapping
    # labels as balanced blocks until they fit. Marker-to-label leader lines
    # preserve identity, so displacement alone must not force rotation.
    my $chosen = $layout eq 'auto'
        ? ($horizontal ? 'horizontal' : 'vertical') : $layout;
    my $centers = $horizontal;
    my $reason = $chosen eq 'horizontal' ? 'horizontal_labels_fit' :
        (!$horizontal ? 'horizontal_width_exceeds_canvas' : 'explicit_vertical');
    if ($chosen eq 'vertical') {
        my @vertical_half = map { $font * 96 / 72 * 0.72 } @items;
        $centers = place_centers(\@items, \@vertical_half, $left, $right, 5);
        die "Too many target labels for the local GTF plot width; increase gtf_design_width\n"
            unless $centers;
        $reason = 'explicit_vertical' if $layout eq 'vertical';
    } elsif (!$centers) {
        die "Horizontal target labels cannot fit; use auto/vertical or increase gtf_design_width\n";
    }
    my @positions;
    for my $i (0 .. $#items) {
        my $bp = int($x_min + ($centers->[$i] - $left) / $plot_width * $span + 0.5);
        push @positions, { snp => $items[$i]{snp}, bp => $bp };
    }
    my $text_height = $chosen eq 'horizontal' ? $font * 96 / 72 * 1.2
        : 2 * (sort { $b <=> $a } @horizontal_half)[0];
    my $headroom_px = $text_height + ($font * 96 / 72 * 0.8 > 12
        ? $font * 96 / 72 * 0.8 : 12);
    # The GTL plot uses the full design height when converting offsetmax to
    # visible headroom. Dividing by a smaller nominal area over-reserves space,
    # especially for rotated labels.
    my $headroom_frac = $headroom_px / $height;
    $headroom_frac = 0.025 if $headroom_frac < 0.025;
    $headroom_frac = 0.45 if $headroom_frac > 0.45;
    my $center_offset = $headroom_frac * $height / (2 * $text_height);
    return { layout => $chosen,
        positions => join(' ', map { "$_->{snp}=$_->{bp}" } @positions),
        font_size => $font, reason => $reason,
        headroom_frac => sprintf('%.4f', $headroom_frac),
        center_offset => sprintf('%.3f', $center_offset),
        max_horizontal_shift_px => int($max_horizontal_shift + 0.5),
        usable_width_px => int($plot_width + 0.5) };
}

unless (caller) {
    my %opt = (layout => 'auto', design_width => 950, font_size => 10);
    GetOptions('targets=s' => \$opt{targets}, 'reference-bp=f' => \$opt{reference_bp},
        'window-bp=f' => \$opt{window_bp}, 'design-width=f' => \$opt{design_width},
        'design-height=f' => \$opt{design_height},
        'font-size=f' => \$opt{font_size}, 'layout=s' => \$opt{layout},
        'chromosome-length=f' => \$opt{chromosome_length}) or die "Invalid label planner option\n";
    die "--targets SNP=BP[,SNP=BP...] is required\n" unless length($opt{targets} // '');
    my @targets = map {
        my ($snp, $bp) = split /=/, $_, 2;
        +{ snp => $snp, bp => $bp }
    } split /,/, $opt{targets};
    print encode_json(plan(%opt, targets => \@targets)), "\n";
}
1;
