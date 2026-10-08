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
ok($plan->{headroom_frac} > 0.03 && $plan->{center_offset} > 0,
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
is($plan->{layout}, 'vertical', 'crowded long labels rotate');
is(scalar(split / /, $plan->{positions}), 8, 'every crowded target receives a position');

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
