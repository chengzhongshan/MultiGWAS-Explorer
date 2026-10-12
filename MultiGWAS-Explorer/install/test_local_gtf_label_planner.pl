#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Test::More;
require "$Bin/../DiffGWASDeps/plan_local_gtf_labels.pl";

my @nearby = (
    { snp => 'rs2070788', bp => 1_000_000 },
    { snp => 'rs383510', bp => 1_000_100 },
);
my $plan = LocalGTFLabelPlanner::plan(targets => \@nearby,
    reference_bp => 1_000_000, window_bp => 1_000_000,
    design_width => 950, font_size => 10);
is($plan->{layout}, 'horizontal', 'two nearby rsIDs fit after balanced displacement');
ok($plan->{headroom_frac} >= 0.025 && $plan->{center_offset} > 0,
    'horizontal labels receive height and centering parameters');
my @adjusted = map { (split /=/, $_, 2)[1] } split / /, $plan->{positions};
ok($adjusted[0] < 1_000_000 && $adjusted[1] > 1_000_100,
    'planner moves labels in both directions');
my $scale = $plan->{usable_width_px} / 2_000_000;
my $distance_px = ($adjusted[1] - $adjusted[0]) * $scale;
my $needed_px = (LocalGTFLabelPlanner::text_width_px('rs2070788', 10)
    + LocalGTFLabelPlanner::text_width_px('rs383510', 10)) / 2
    + 10 * 96 / 72 * 0.7 + 6;
ok($distance_px + 0.1 >= $needed_px,
    'adjusted horizontal labels have a font-aware gap');

my @crowded = map { +{ snp => "rs123456789$_", bp => 1_000_000 + $_ } } 1 .. 8;
$plan = LocalGTFLabelPlanner::plan(targets => \@crowded,
    reference_bp => 1_000_000, window_bp => 1_000_000,
    design_width => 950, font_size => 10);
is($plan->{layout}, 'vertical', 'labels rotate when their total width cannot fit');
is(scalar(split / /, $plan->{positions}), 8, 'every crowded target receives a position');

my @cd55_targets = (
    { snp => 'rs2802216', bp => 207_394_658 },
    { snp => 'rs200078770', bp => 207_394_677 },
    { snp => 'rs4303075', bp => 207_394_678 },
    { snp => 'rs11117664', bp => 207_394_689 },
);
$plan = LocalGTFLabelPlanner::plan(targets => \@cd55_targets,
    reference_bp => 207_394_658, window_bp => 1_000_031,
    design_width => 950, design_height => 1000, font_size => 10);
is($plan->{layout}, 'horizontal', 'four CD55-region targets fit after block spacing');
is(scalar(split / /, $plan->{positions}), 4,
    'four CD55-region targets all retain labels');
ok($plan->{headroom_frac} < 0.05,
    'four adjusted horizontal labels use compact headroom');

my @chr21_targets = (
    { snp => 'rs2070788', bp => 41_470_061 },
    { snp => 'rs383510', bp => 41_486_440 },
    { snp => 'rs462687', bp => 41_425_130 },
    { snp => 'rs429442', bp => 41_489_405 },
);
$plan = LocalGTFLabelPlanner::plan(targets => \@chr21_targets,
    reference_bp => 41_470_061, window_bp => 1_044_931,
    design_width => 950, design_height => 1000, font_size => 10);
is($plan->{layout}, 'horizontal', 'four chr21 targets fit horizontally with balanced leader lines');
is(scalar(split / /, $plan->{positions}), 4, 'every chr21 target retains a label');
ok($plan->{max_horizontal_shift_px} > 90,
    'a large but feasible chr21 shift does not force rotation');
my @chr21_placed = map { [split /=/, $_, 2] } split / /, $plan->{positions};
my $chr21_px_per_bp = $plan->{usable_width_px} / (2 * 1_044_931);
my $chr21_gap_px = 10 * 96 / 72 * 0.7 + 6;
for my $i (1 .. $#chr21_placed) {
    my ($left_snp, $left_bp) = @{ $chr21_placed[$i - 1] };
    my ($right_snp, $right_bp) = @{ $chr21_placed[$i] };
    my $required_px = (LocalGTFLabelPlanner::text_width_px($left_snp, 10)
        + LocalGTFLabelPlanner::text_width_px($right_snp, 10)) / 2 + $chr21_gap_px;
    ok(($right_bp - $left_bp) * $chr21_px_per_bp + 0.1 >= $required_px,
        "adjusted labels $left_snp and $right_snp do not overlap");
}

my @medium = map { +{ snp => "rs12345$_", bp => 1_000_000 + $_ * 120_000 } } 1 .. 5;
my $small = LocalGTFLabelPlanner::plan(targets => \@medium,
    reference_bp => 1_000_000, window_bp => 1_000_000,
    design_width => 950, font_size => 7);
my $large = LocalGTFLabelPlanner::plan(targets => \@medium,
    reference_bp => 1_000_000, window_bp => 1_000_000,
    design_width => 950, font_size => 20);
is($small->{layout}, 'horizontal', 'smaller font fits the same target set');
is($large->{layout}, 'vertical', 'larger font changes the decision');
ok($large->{headroom_frac} > $small->{headroom_frac},
    'larger vertical text receives more top headroom');

$plan = LocalGTFLabelPlanner::plan(targets => \@nearby,
    reference_bp => 1_000_000, window_bp => 1_000_000,
    design_width => 950, font_size => 10, layout => 'vertical');
is($plan->{layout}, 'vertical', 'explicit vertical layout is respected');

my $failed = eval {
    LocalGTFLabelPlanner::plan(targets => \@crowded,
        reference_bp => 1_000_000, window_bp => 1_000_000,
        design_width => 950, font_size => 10, layout => 'horizontal');
    1;
};
ok(!$failed && $@ =~ /Horizontal target labels cannot fit/,
    'impossible forced horizontal layout gives an actionable error');

done_testing();
