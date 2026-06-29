#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $input =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.tsv.gz';
my $output =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.bgz.tsv.gz';
my $htsbin =
  '/mnt/g/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/local/bin';
my $seq_col   = 1;
my $start_col = 2;
my $end_col   = 2;

GetOptions(
    'input=s'  => \$input,
    'output=s' => \$output,
    'htsbin=s' => \$htsbin,
    'seq=i'    => \$seq_col,
    'start=i'  => \$start_col,
    'end=i'    => \$end_col,
) or die usage();

my $bgzip = -x "$htsbin/bgzip" ? "$htsbin/bgzip" : 'bgzip';
my $tabix = -x "$htsbin/tabix" ? "$htsbin/tabix" : 'tabix';

die "Input file not found: $input\n" unless -s $input;

open my $in, '-|', "zcat '$input'" or die "Cannot read $input with zcat: $!\n";
open my $out, '|-', "'$bgzip' -@ 4 -c > '$output'"
  or die "Cannot write bgzip output $output: $!\n";

my $header = <$in>;
die "Input is empty: $input\n" unless defined $header;
chomp $header;
$header =~ s/\r$//;
$header =~ s/^#//;
print {$out} "#$header\n";

my $rows = 0;
while (my $line = <$in>) {
    print {$out} $line;
    $rows++;
}

close $in;
close $out or die "Failed closing bgzip output $output: $!\n";

system($tabix, '-f', '-s', $seq_col, '-b', $start_col, '-e', $end_col, '-S', 1, $output) == 0
  or die "tabix failed for $output\n";

print "Input:  $input\n";
print "Output: $output\n";
print "Index:  $output.tbi\n";
print "Rows:   $rows\n";

sub usage {
    return <<"USAGE";
Usage:
  perl bgzip_tabix_diff_gwas.pl [options]

Options:
  --input FILE.tsv.gz       Sorted gzip input table
  --output FILE.tsv.gz      bgzip output table
  --htsbin DIR              Directory containing bgzip/tabix
  --seq N                   1-based chromosome column. Default: 1
  --start N                 1-based start column. Default: 2
  --end N                   1-based end column. Default: 2
USAGE
}
