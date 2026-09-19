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
my $im=GD::Image->new($fh) or die 'Invalid PNG'; close $fh;
# Exclude the plot borders and data area: decoding alone cannot detect lost text.
cmp_ok(ink($im,20,40,140,360),'>',20,'cohort labels are visible in the left margin');
cmp_ok(ink($im,220,385,680,410),'>',20,'x-axis title is visible below the plot');
open $csv,'>',"$dir/multi.csv" or die $!;
print {$csv} "SNP,gene,hit_class,BETA,SE,P\n",
 "rs1,ZSCAN12,DIFFERENTIAL,0.1,0.05,0.04\n",
 "rs2,ZSCAN12,DIFFERENTIAL,0.2,0.05,0.01\n",
 "rs3,TENM1,DIFFERENTIAL,-0.1,0.05,0.04\n";
close $csv;
$rc=system($^X,"$Bin/../DiffGWASDeps/gnuplot/pdl_gunplot_forest.pl",
 '--csv',"$dir/multi.csv",'--out-prefix',"$dir/multi",'--track-ids','TEST',
 '--track-labels','Test_cohort','--track-beta-vars','BETA','--track-se-vars','SE',
 '--track-p-vars','P','--width',900,'--height',420);
is($rc,0,'multi-SNP forest renderer succeeds');
open my $gp,'<',"$dir/multi_TEST.gp" or die $!;
my $script=do {local $/; <$gp>}; close $gp;
like($script,qr/set y2tics[^\n]*"ZSCAN12" 1, "ZSCAN12" 2, "TENM1" 3/,
 'right-axis gene symbols retain distinct SNP row coordinates');
open $fh,'<:raw',"$dir/multi_TEST.png" or die $!;
$im=GD::Image->new($fh) or die 'Invalid multi-SNP PNG'; close $fh;
for my $band ([40,140],[140,250],[250,360]) {
 cmp_ok(ink($im,730,$band->[0],890,$band->[1]),'>',20,
  "gene label visible in right-margin row $band->[0]-$band->[1]");
}
done_testing;
sub ink {
 my ($im,$x0,$y0,$x1,$y1)=@_; my $n=0;
 for my $x($x0..$x1) {for my $y($y0..$y1) {
  my @rgb=$im->rgb($im->getPixel($x,$y)); ++$n if $rgb[0]<180 && $rgb[1]<180 && $rgb[2]<180;
 }} return $n;
}
