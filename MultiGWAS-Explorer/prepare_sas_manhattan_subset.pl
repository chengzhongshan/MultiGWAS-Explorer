#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $input =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.tsv.gz';
my $output =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.SAS_manhattan_p_lt_0p05.tsv.gz';
my $manifest =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.SAS_manhattan_p_lt_0p05.manifest.tsv';
my $threshold = 0.05;

GetOptions(
    'input=s'     => \$input,
    'output=s'    => \$output,
    'manifest=s'  => \$manifest,
    'threshold=f' => \$threshold,
) or die usage();

die "Input file not found: $input\n" unless -s $input;
die "threshold must be positive\n" unless $threshold > 0;

my @out_cols = qw(
  CHR BP SNP
  ALL_STD_P ASN_STD_P EUR_STD_P
  ALL_DIFF_P ASN_DIFF_P EUR_DIFF_P
);

my %pair_to_prefix = (
    SCZ_W3_ALL_SEX => 'ALL',
    SCZ_W3_ASN_SEX => 'ASN',
    SCZ_W3_EUR_SEX => 'EUR',
);

open my $in, '-|', "zcat '$input'" or die "Cannot read $input with zcat: $!\n";
open my $out, '|-', "gzip -c > '$output'" or die "Cannot write $output with gzip: $!\n";
print {$out} join("\t", @out_cols), "\n";

my $header = <$in>;
die "Input is empty: $input\n" unless defined $header;
chomp $header;
$header =~ s/\r$//;
$header =~ s/^#//;
my @cols = split /\t/, $header, -1;
my %idx;
for my $i (0 .. $#cols) {
    $idx{$cols[$i]} = $i;
}
for my $required (qw(CHR BP SNP PAIR_TAG DIFF_P STD_DIFF_P)) {
    die "Required column $required not found in header\n" unless exists $idx{$required};
}

my %stats = (
    rows_read      => 0,
    groups_seen    => 0,
    rows_written   => 0,
    groups_skipped => 0,
    bad_p          => 0,
);

my $current_key = '';
my @bucket;
while (my $line = <$in>) {
    chomp $line;
    $line =~ s/\r$//;
    next if $line =~ /^\s*$/;
    my @v = split /\t/, $line, -1;
    my $key = join("\t", map { $v[ $idx{$_} ] // '' } qw(CHR BP SNP));
    if (@bucket && $key ne $current_key) {
        process_bucket(\@bucket, $out);
        @bucket = ();
    }
    $current_key = $key;
    push @bucket, \@v;
    $stats{rows_read}++;
}
process_bucket(\@bucket, $out) if @bucket;

close $in;
close $out or die "Failed closing gzip output $output: $!\n";

open my $man, '>', $manifest or die "Cannot write $manifest: $!\n";
print {$man} join("\t", qw(METRIC VALUE)), "\n";
for my $metric (sort keys %stats) {
    print {$man} join("\t", $metric, $stats{$metric}), "\n";
}
print {$man} join("\t", 'input', $input), "\n";
print {$man} join("\t", 'output', $output), "\n";
print {$man} join("\t", 'threshold', $threshold), "\n";
print {$man} join("\t", 'columns', join(',', @out_cols)), "\n";
close $man;

print "Input:        $input\n";
print "Output:       $output\n";
print "Manifest:     $manifest\n";
print "Rows read:    $stats{rows_read}\n";
print "Rows written: $stats{rows_written}\n";

sub process_bucket {
    my ($bucket, $out) = @_;
    $stats{groups_seen}++;

    my %row = map { $_ => '' } @out_cols;
    $row{CHR} = sas_chr($bucket->[0][ $idx{CHR} ] // '');
    $row{BP}  = $bucket->[0][ $idx{BP} ]  // '';
    $row{SNP} = $bucket->[0][ $idx{SNP} ] // '';

    my $keep = 0;
    for my $v (@$bucket) {
        my $pair = $v->[ $idx{PAIR_TAG} ] // '';
        my $prefix = $pair_to_prefix{$pair} // next;

        my $std_p  = numeric($v->[ $idx{STD_DIFF_P} ]);
        my $diff_p = numeric($v->[ $idx{DIFF_P} ]);

        if (defined $std_p) {
            $row{"${prefix}_STD_P"} = $std_p;
            $keep = 1 if $std_p < $threshold;
        }
        else {
            $stats{bad_p}++;
        }

        if (defined $diff_p) {
            $row{"${prefix}_DIFF_P"} = $diff_p;
            $keep = 1 if $diff_p < $threshold;
        }
        else {
            $stats{bad_p}++;
        }
    }

    if ($keep) {
        print {$out} join("\t", map { defined $row{$_} ? $row{$_} : '' } @out_cols), "\n";
        $stats{rows_written}++;
    }
    else {
        $stats{groups_skipped}++;
    }
}

sub numeric {
    my ($x) = @_;
    return undef unless defined $x;
    return undef if $x eq '' || $x =~ /^(?:NA|NaN|null|\.)$/i;
    return undef unless $x =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
    return 0 + $x;
}

sub sas_chr {
    my ($chr) = @_;
    $chr =~ s/^chr//i;
    return 23 if $chr =~ /^(?:X|23)$/i;
    return $chr;
}

sub usage {
    return <<"USAGE";
Usage:
  perl prepare_sas_manhattan_subset.pl [options]

Options:
  --input FILE.tsv.gz       Standardized differential GWAS table
  --output FILE.tsv.gz      Small wide-format SAS Manhattan input
  --manifest FILE.tsv       Run summary
  --threshold FLOAT         Keep SNP if any selected P column < threshold. Default: 0.05

Output columns:
  CHR BP SNP ALL_STD_P ASN_STD_P EUR_STD_P ALL_DIFF_P ASN_DIFF_P EUR_DIFF_P
USAGE
}
