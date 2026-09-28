#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use File::Spec;
use lib "$Bin/DiffGWASDeps";
use HTSToolResolver qw(resolve_hts_tool external_tool_path);

my $input =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.tsv.gz';
my $output =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.bgz.tsv.gz';
my $htsbin = File::Spec->catdir($Bin, 'local', 'bin');
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

my $bgzip = resolve_hts_tool('bgzip', explicit => $htsbin, start_dir => $Bin);
my $tabix = resolve_hts_tool('tabix', explicit => $htsbin, start_dir => $Bin);

die "Native bgzip/tabix were not found. On Cygwin run bash install/repair_and_test_cygwin.sh; "
  . "Windows executables inherited from the global PATH are intentionally ignored.\n"
  unless $bgzip && $tabix;

die "Input file not found: $input\n" unless -s $input;

open my $in, '-|', 'gzip', '-dc', '--', $input or die "Cannot read $input with gzip: $!\n";
pipe(my $bgzip_reader, my $out) or die "Cannot create bgzip pipe: $!\n";
my $bgzip_pid = fork();
die "Cannot fork bgzip: $!\n" unless defined $bgzip_pid;
if ($bgzip_pid == 0) {
  close $out;
  open STDIN, '<&', $bgzip_reader or die "Cannot connect bgzip stdin: $!\n";
  open STDOUT, '>', $output or die "Cannot open bgzip output $output: $!\n";
  exec { $bgzip } $bgzip, '-c';
  die "Cannot execute $bgzip: $!\n";
}
close $bgzip_reader;

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
close $out or die "Failed closing bgzip input for $output: $!\n";
waitpid($bgzip_pid, 0);
die "bgzip failed for $output\n" unless $? == 0 && -s $output;

system($tabix, '-f', '-s', $seq_col, '-b', $start_col, '-e', $end_col, '-S', 1,
  external_tool_path($tabix, $output)) == 0
  or die "tabix failed for $output\n";

die "tabix did not create the expected index: $output.tbi\n" unless -s "$output.tbi";

print "Input:  $input\n";
print "Output: $output\n";
print "Index:  $output.tbi\n";
print "Rows:   $rows\n";
print "bgzip:  $bgzip\n";
print "tabix:  $tabix\n";

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
