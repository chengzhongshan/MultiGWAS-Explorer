#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use File::Spec;

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

my $bgzip = resolve_hts_tool('bgzip', $htsbin);
my $tabix = resolve_hts_tool('tabix', $htsbin);

die "Input file not found: $input\n" unless -s $input;

open my $in, '-|', "zcat '$input'" or die "Cannot read $input with zcat: $!\n";
my $bgzip_args = $bgzip =~ /\.exe$/i ? '-c' : '-@ 4 -c';
open my $out, '|-', "'$bgzip' $bgzip_args > '$output'"
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

my $tabix_output_arg = path_for_tool($tabix, $output);
system($tabix, '-f', '-s', $seq_col, '-b', $start_col, '-e', $end_col, '-S', 1, $tabix_output_arg) == 0
  or die "tabix failed for $output\n";

die "tabix did not create the expected index: $output.tbi\n" unless -s "$output.tbi";

print "Input:  $input\n";
print "Output: $output\n";
print "Index:  $output.tbi\n";
print "Rows:   $rows\n";
print "bgzip:  $bgzip\n";
print "tabix:  $tabix\n";

sub resolve_hts_tool {
    my ($name, $requested_dir) = @_;
    my @dirs = grep { defined($_) && length($_) } (
        $requested_dir,
        File::Spec->catdir($Bin, 'DiffGWASDeps'),
    );
    my %seen;
    for my $dir (@dirs) {
        next if $seen{$dir}++;
        my @tool_names = $^O eq 'cygwin' ? ("$name.exe", $name) : ($name, "$name.exe");
        for my $tool_name (@tool_names) {
            my $candidate = File::Spec->catfile($dir, $tool_name);
            return $candidate if -f $candidate
              && (-x $candidate || ($^O eq 'cygwin' && $candidate =~ /\.exe$/i));
        }
    }
    return $name;
}

sub path_for_tool {
    my ($tool, $path) = @_;
    return $path unless $^O eq 'cygwin' && $tool =~ /\.exe$/i;

    open my $cygpath, '-|', 'cygpath', '-w', $path
      or die "Cannot convert Cygwin path for Windows executable $tool: $!\n";
    my $windows_path = <$cygpath>;
    close $cygpath
      or die "cygpath failed while converting $path for $tool\n";
    die "cygpath returned no Windows path for $path\n"
      unless defined($windows_path) && length($windows_path);
    $windows_path =~ s/[\r\n]+\z//;
    return $windows_path;
}

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
