#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);

my $input =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.tsv.gz';
my $output =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.tsv.gz';
my $manifest =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.manifest.tsv';
my $z_col = 'DIFF_Z';
my $htsbin = "$Bin/../local/bin";
my $index_output = 1;

GetOptions(
    'input=s'    => \$input,
    'output=s'   => \$output,
    'manifest=s' => \$manifest,
    'z-col=s'    => \$z_col,
    'htsbin=s'   => \$htsbin,
    'index-output!' => \$index_output,
) or die usage();

die "Input file not found: $input\n" unless -s $input;

my ($header, $cols, $idx) = read_header($input);
my %idx = %$idx;
die "Column $z_col not found in $input\n" unless exists $idx{$z_col};
my $bgzip = resolve_hts_tool($htsbin, 'bgzip');
my $tabix = resolve_hts_tool($htsbin, 'tabix');
my $can_index_output = (
    $index_output
    && defined $bgzip
    && defined $tabix
    && $output =~ /\.gz$/i
    && exists $idx{CHR}
    && exists $idx{BP}
);

my ($n, $mean, $sd) = z_stats($input, $idx{$z_col});
die "No numeric values found in $z_col\n" unless $n > 0;
die "Cannot standardize because standard deviation is zero\n" unless $sd > 0;

open my $in,  '-|', "zcat '$input'"       or die "Cannot read $input with zcat: $!\n";
my $writer_cmd = $can_index_output
  ? "$bgzip -@ 4 -c > '$output'"
  : "gzip -c > '$output'";
open my $out, '|-', $writer_cmd or die "Cannot write $output: $!\n";

my $raw_header = <$in>;
chomp $raw_header;
$raw_header =~ s/\r$//;
print {$out} join("\t", $raw_header, 'ORIG_DIFF_Z', 'ORIG_DIFF_P', 'STD_DIFF_Z', 'STD_DIFF_P'), "\n";

my $rows = 0;
my $numeric_rows = 0;
while (my $line = <$in>) {
    chomp $line;
    $line =~ s/\r$//;
    my @v = split /\t/, $line, -1;
    my $z = numeric($v[ $idx{$z_col} ]);
    if (defined $z) {
        my $std_z = ($z - $mean) / $sd;
        my $std_p = two_sided_p_from_z($std_z);
        print {$out} join("\t", @v, fmt($z), value(\@v, \%idx, 'DIFF_P'), fmt($std_z), p_fmt($std_p)), "\n";
        $numeric_rows++;
    }
    else {
        print {$out} join("\t", @v, '', value(\@v, \%idx, 'DIFF_P'), '', ''), "\n";
    }
    $rows++;
}

close $in;
close $out or die "Failed closing gzip output $output: $!\n";

my $index_status = 'disabled';
if ($can_index_output) {
    my $seq_col = $idx{CHR} + 1;
    my $bp_col  = $idx{BP} + 1;
    if (system($tabix, '-f', '-s', $seq_col, '-b', $bp_col, '-e', $bp_col, '-S', 1, $output) == 0) {
        $index_status = 'created';
    }
    else {
        warn "Warning: tabix failed for $output; continuing without an index.\n";
        unlink "$output.tbi" if -e "$output.tbi";
        unlink "$output.csi" if -e "$output.csi";
        $index_status = 'failed';
    }
}
elsif ($index_output) {
    $index_status = 'unavailable';
}

open my $man, '>', $manifest or die "Cannot write $manifest: $!\n";
print {$man} join("\t", qw(METRIC VALUE)), "\n";
print {$man} join("\t", 'input',        $input), "\n";
print {$man} join("\t", 'output',       $output), "\n";
print {$man} join("\t", 'z_col',        $z_col), "\n";
print {$man} join("\t", 'index_output', $index_output ? 1 : 0), "\n";
print {$man} join("\t", 'index_status', $index_status), "\n";
print {$man} join("\t", 'numeric_n',    $n), "\n";
print {$man} join("\t", 'mean',         fmt($mean)), "\n";
print {$man} join("\t", 'sample_sd',    fmt($sd)), "\n";
print {$man} join("\t", 'rows_written', $rows), "\n";
print {$man} join("\t", 'numeric_rows', $numeric_rows), "\n";
close $man;

print "Input:        $input\n";
print "Output:       $output\n";
print "Manifest:     $manifest\n";
print "Index status: $index_status\n";
print "Numeric N:    $n\n";
print "Mean $z_col:  ", fmt($mean), "\n";
print "SD $z_col:    ", fmt($sd), "\n";
print "Rows written: $rows\n";
print "Index:        $output.tbi\n" if $index_status eq 'created';

sub resolve_hts_tool {
    my ($dir, $tool) = @_;
    for my $candidate (
        (defined $dir && length $dir ? ("$dir/$tool", "$dir/$tool.exe") : ()),
        $tool,
    ) {
        next unless defined $candidate && length $candidate;
        return $candidate if -x $candidate || command_exists($candidate);
    }
    return undef;
}

sub command_exists {
    my ($cmd) = @_;
    return 0 unless defined $cmd && length $cmd;
    return scalar(`command -v '$cmd' 2>/dev/null`) ? 1 : 0;
}

sub read_header {
    my ($path) = @_;
    open my $fh, '-|', "zcat '$path'" or die "Cannot read $path with zcat: $!\n";
    my $h = <$fh>;
    close $fh;
    die "Input is empty: $path\n" unless defined $h;
    chomp $h;
    $h =~ s/\r$//;
    $h =~ s/^#//;
    my @c = split /\t/, $h, -1;
    my %i;
    for my $j (0 .. $#c) {
        $i{$c[$j]} = $j;
    }
    return ($h, \@c, \%i);
}

sub z_stats {
    my ($path, $col_i) = @_;
    open my $fh, '-|', "zcat '$path'" or die "Cannot read $path with zcat: $!\n";
    <$fh>;
    my ($n, $mean, $m2) = (0, 0, 0);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        my @v = split /\t/, $line, -1;
        my $x = numeric($v[$col_i]);
        next unless defined $x;
        $n++;
        my $delta = $x - $mean;
        $mean += $delta / $n;
        my $delta2 = $x - $mean;
        $m2 += $delta * $delta2;
    }
    close $fh;
    my $sd = $n > 1 ? sqrt($m2 / ($n - 1)) : 0;
    return ($n, $mean, $sd);
}

sub value {
    my ($vals, $idx, $col) = @_;
    return '' unless exists $idx->{$col};
    return $vals->[ $idx->{$col} ] // '';
}

sub numeric {
    my ($x) = @_;
    return undef unless defined $x;
    return undef if $x eq '' || $x =~ /^(?:NA|NaN|null|\.)$/i;
    return undef unless $x =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
    return 0 + $x;
}

sub two_sided_p_from_z {
    my ($z) = @_;
    return erfc_approx(abs($z) / sqrt(2));
}

sub erfc_approx {
    my ($x) = @_;
    my $z = abs($x);
    my $t = 1 / (1 + 0.5 * $z);
    my $poly = ((((((((0.17087277 * $t - 0.82215223) * $t + 1.48851587) * $t
      - 1.13520398) * $t + 0.27886807) * $t - 0.18628806) * $t
      + 0.09678418) * $t + 0.37409196) * $t + 1.00002368) * $t;
    my $r = $t * exp(-$z * $z - 1.26551223 + $poly);
    return $x >= 0 ? $r : 2 - $r;
}

sub fmt {
    my ($x) = @_;
    return '' unless defined $x;
    return sprintf('%.10g', $x);
}

sub p_fmt {
    my ($x) = @_;
    return '' unless defined $x;
    return sprintf('%.6e', $x);
}

sub usage {
    return <<"USAGE";
Usage:
  perl standardize_diff_gwas_zscore.pl [options]

Options:
  --input FILE.tsv.gz       Differential GWAS table
  --output FILE.tsv.gz      Output with ORIG_DIFF_Z, ORIG_DIFF_P, STD_DIFF_Z, STD_DIFF_P
  --manifest FILE.tsv       Run summary with mean and sample SD
  --z-col NAME              Z-score column to standardize. Default: DIFF_Z
  --htsbin DIR              Directory containing bgzip/tabix. Default: ../local/bin
  --index-output            Try to bgzip-index the standardized output with tabix. Default: on
  --no-index-output         Write plain gzip output without tabix indexing
USAGE
}
