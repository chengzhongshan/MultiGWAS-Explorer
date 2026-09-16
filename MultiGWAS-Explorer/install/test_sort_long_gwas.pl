#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);
my $dir=tempdir('multigwas sort XXXXX',TMPDIR=>1,CLEANUP=>1);
my $input="CHR\tBP\tSNP\n23\t200\trs3\n1\t20\trs2\n1\t10\trs1\n1\tNA\trs4\n";
gzip(\$input=>"$dir/input.gz") or die $GzipError;
local $ENV{INPUT_GZ}="$dir/input.gz";
local $ENV{OUTPUT_GZ}="$dir/sorted.gz";
local $ENV{EXCLUDED_GZ}="$dir/excluded.gz";
local $ENV{TMPDIR_SORT}="$dir/tmp";
system('bash',"$Bin/../DiffGWASDeps/sort_long_gwas_by_coord.sh")==0 or die "Sorter failed\n";
my $fh=IO::Uncompress::Gunzip->new("$dir/sorted.gz",MultiStream=>1) or die $GunzipError;
my @lines=<$fh>;
$fh->close or die "Cannot close sorted gzip\n";
die "Wrong sorted content\n" unless join('',@lines) eq "#CHR\tBP\tSNP\n1\t10\trs1\n1\t20\trs2\n23\t200\trs3\n";
my $indexed="$dir/sorted.gz";
if ($^O eq 'cygwin') {
 open my $cp,'-|','cygpath','-w',$indexed or die $!;
 $indexed=<$cp>; chomp $indexed; close $cp or die 'cygpath failed';
}
for my $region ('1:1-30','23:1-300') {
 open my $tab,'-|','tabix',$indexed,$region or die $!;
 my @rows=<$tab>;
 close $tab or die "tabix query failed\n";
 die "Empty region $region\n" unless @rows;
 die "Wrong query rows for $region\n" unless @rows==($region=~/^1:/?2:1);
}
print "PASS: coordinate sorting, excluded row, real tabix index, autosome/X queries, paths containing spaces\n";
my $diff="CHR\tBP\tSNP\tDIFF_Z\tDIFF_P\n1\t10\trs1\t1\t0.3\n1\t20\trs2\t-1\t0.3\n23\t200\trs3\t2\t0.05\n";
gzip(\$diff=>"$dir/diff.gz") or die $GzipError;
system($^X,"$Bin/../DiffGWASDeps/standardize_diff_gwas_zscore.pl",
 '--input',"$dir/diff.gz",'--output',"$dir/std.gz",'--manifest',"$dir/std.manifest.tsv")==0
 or die "Standardization failed\n";
open my $mf,'<',"$dir/std.manifest.tsv" or die $!;
my $metrics=do {local $/;<$mf>};close $mf;
die "Standardized output was not indexed\n" unless $metrics=~/^index_status\tcreated$/m;
print "PASS: standardized differential output has a tabix index\n";
