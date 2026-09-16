#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use GD;
use Test::More;
my $dir=tempdir('forest-text-XXXXX',TMPDIR=>1,CLEANUP=>1);
open my $csv,'>',"$dir/hits.csv" or die $!;
print {$csv} "SNP,hit_class,BETA,SE,P\nrs_123*,DIFFERENTIAL,0.1,0.05,0.04\n";
close $csv;
my $rc=system($^X,"$Bin/../DiffGWASDeps/gnuplot/pdl_gunplot_forest.pl",
 '--csv',"$dir/hits.csv",'--out-prefix',"$dir/forest",'--track-ids','TEST',
 '--track-labels','Test_cohort','--track-beta-vars','BETA','--track-se-vars','SE',
 '--track-p-vars','P','--width',900,'--height',420);
is($rc,0,'forest renderer succeeds');
open my $fh,'<:raw',"$dir/forest_single_snp.png" or die $!;
my $im=GD::Image->newFromPng($fh,1) or die 'Invalid PNG'; close $fh;
# Exclude the plot borders and data area: decoding alone cannot detect lost text.
cmp_ok(ink($im,20,40,140,360),'>',20,'cohort labels are visible in the left margin');
cmp_ok(ink($im,220,385,680,410),'>',20,'x-axis title is visible below the plot');
done_testing;
sub ink {
 my ($im,$x0,$y0,$x1,$y1)=@_; my $n=0;
 for my $x($x0..$x1) {for my $y($y0..$y1) {
  my @rgb=$im->rgb($im->getPixel($x,$y)); ++$n if $rgb[0]<180 && $rgb[1]<180 && $rgb[2]<180;
 }} return $n;
}
