#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $input =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_sex_stratified_merged_long.sorted.coord.tsv.gz';
my $output =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.tsv.gz';
my $manifest =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.manifest.tsv';
my $beta_col = '';
my $or_col   = 'OR';
my $se_col   = 'SE';
my $rho      = 0;
my $limit    = 0;

GetOptions(
    'input=s'    => \$input,
    'output=s'   => \$output,
    'manifest=s' => \$manifest,
    'beta-col=s' => \$beta_col,
    'or-col=s'   => \$or_col,
    'se-col=s'   => \$se_col,
    'rho=f'      => \$rho,
    'limit=i'    => \$limit,
) or die usage();

die usage() if !defined $input || !defined $output || !defined $manifest;
die "rho must be greater than -1 and less than 1\n" unless $rho > -1 && $rho < 1;

my @out_cols = qw(
  CHR BP A1 A2 SNP PAIR_TAG
  FEMALE_GWAS_TAG MALE_GWAS_TAG
  FEMALE_SOURCE_FILE MALE_SOURCE_FILE
  FEMALE_BETA MALE_BETA DIFF_BETA
  FEMALE_SE MALE_SE DIFF_SE
  FEMALE_Z MALE_Z DIFF_Z DIFF_P
  FEMALE_OR MALE_OR FEMALE_P MALE_P
  FEMALE_FRQ_A FEMALE_FRQ_U MALE_FRQ_A MALE_FRQ_U
  FEMALE_INFO MALE_INFO
  CHR_ORIGINAL IS_CHRX
);

my %stats = (
    rows_read       => 0,
    groups_seen     => 0,
    pairs_written   => 0,
    skipped_no_pair => 0,
    skipped_bad_num => 0,
    skipped_bad_var => 0,
);

open my $in,  '-|', "zcat '$input'"       or die "Cannot read $input with zcat: $!\n";
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

for my $required (qw(CHR BP A1 A2 SNP GWAS_TAG SOURCE_FILE CHR_ORIGINAL IS_CHRX SE P)) {
    die "Required column $required not found in header: $header\n" unless exists $idx{$required};
}
if ($beta_col) {
    die "Requested beta column $beta_col not found in header\n" unless exists $idx{$beta_col};
}
else {
    die "OR column $or_col not found in header, and --beta-col was not provided\n"
      unless exists $idx{$or_col};
}

my $current_key = '';
my @bucket;

while (my $line = <$in>) {
    chomp $line;
    $line =~ s/\r$//;
    next if $line =~ /^\s*$/;

    my @v = split /\t/, $line, -1;
    my $key = join("\t", map { value(\@v, $_) } qw(CHR BP A1 A2 SNP));

    if (@bucket && $key ne $current_key) {
        process_bucket(\@bucket);
        @bucket = ();
        last if $limit && $stats{pairs_written} >= $limit;
    }

    $current_key = $key;
    push @bucket, \@v;
    $stats{rows_read}++;
}
process_bucket(\@bucket) if @bucket && (!$limit || $stats{pairs_written} < $limit);

close $in;
close $out or die "Failed closing gzip output $output: $!\n";

open my $man, '>', $manifest or die "Cannot write $manifest: $!\n";
print {$man} join("\t", qw(METRIC VALUE)), "\n";
for my $metric (sort keys %stats) {
    print {$man} join("\t", $metric, $stats{$metric}), "\n";
}
print {$man} join("\t", 'input',    $input), "\n";
print {$man} join("\t", 'output',   $output), "\n";
print {$man} join("\t", 'beta_col', $beta_col ? $beta_col : "log($or_col)"), "\n";
print {$man} join("\t", 'se_col',   $se_col), "\n";
print {$man} join("\t", 'rho',      $rho), "\n";
close $man;

warn "Differential female-vs-male GWAS output: $output\n";
warn "Manifest: $manifest\n";
warn "Pairs written: $stats{pairs_written}\n";

sub process_bucket {
    my ($bucket) = @_;
    $stats{groups_seen}++;

    my (%female, %male);
    for my $row (@$bucket) {
        my $tag = value($row, 'GWAS_TAG');
        my $pair = pair_tag($tag);
        next unless $pair;

        if ($tag =~ /FEMALE/i) {
            $female{$pair} ||= $row;
        }
        elsif ($tag =~ /MALE/i) {
            $male{$pair} ||= $row;
        }
    }

    my %pairs = map { $_ => 1 } (keys %female, keys %male);
    for my $pair (sort keys %pairs) {
        my $f = $female{$pair};
        my $m = $male{$pair};
        if (!$f || !$m) {
            $stats{skipped_no_pair}++;
            next;
        }

        my $fb = beta($f);
        my $mb = beta($m);
        my $fs = numeric(value($f, $se_col));
        my $ms = numeric(value($m, $se_col));
        if (!defined $fb || !defined $mb || !defined $fs || !defined $ms || $fs <= 0 || $ms <= 0) {
            $stats{skipped_bad_num}++;
            next;
        }

        my $diff_beta = $fb - $mb;
        my $var = $fs * $fs + $ms * $ms - 2 * $rho * $fs * $ms;
        if ($var <= 0) {
            $stats{skipped_bad_var}++;
            next;
        }

        my $diff_se = sqrt($var);
        my $diff_z  = $diff_beta / $diff_se;
        my $diff_p  = two_sided_p_from_z($diff_z);

        my %o = (
            CHR                => value($f, 'CHR'),
            BP                 => value($f, 'BP'),
            A1                 => value($f, 'A1'),
            A2                 => value($f, 'A2'),
            SNP                => value($f, 'SNP'),
            PAIR_TAG           => $pair,
            FEMALE_GWAS_TAG    => value($f, 'GWAS_TAG'),
            MALE_GWAS_TAG      => value($m, 'GWAS_TAG'),
            FEMALE_SOURCE_FILE => value($f, 'SOURCE_FILE'),
            MALE_SOURCE_FILE   => value($m, 'SOURCE_FILE'),
            FEMALE_BETA        => fmt($fb),
            MALE_BETA          => fmt($mb),
            DIFF_BETA          => fmt($diff_beta),
            FEMALE_SE          => value($f, $se_col),
            MALE_SE            => value($m, $se_col),
            DIFF_SE            => fmt($diff_se),
            FEMALE_Z           => fmt($fb / $fs),
            MALE_Z             => fmt($mb / $ms),
            DIFF_Z             => fmt($diff_z),
            DIFF_P             => p_fmt($diff_p),
            FEMALE_OR          => exists $idx{$or_col} ? value($f, $or_col) : '',
            MALE_OR            => exists $idx{$or_col} ? value($m, $or_col) : '',
            FEMALE_P           => value($f, 'P'),
            MALE_P             => value($m, 'P'),
            FEMALE_FRQ_A       => value($f, 'FRQ_A'),
            FEMALE_FRQ_U       => value($f, 'FRQ_U'),
            MALE_FRQ_A         => value($m, 'FRQ_A'),
            MALE_FRQ_U         => value($m, 'FRQ_U'),
            FEMALE_INFO        => value($f, 'INFO'),
            MALE_INFO          => value($m, 'INFO'),
            CHR_ORIGINAL       => value($f, 'CHR_ORIGINAL'),
            IS_CHRX            => value($f, 'IS_CHRX'),
        );

        print {$out} join("\t", map { defined $o{$_} ? $o{$_} : '' } @out_cols), "\n";
        $stats{pairs_written}++;
        last if $limit && $stats{pairs_written} >= $limit;
    }
}

sub pair_tag {
    my ($tag) = @_;
    return '' unless defined $tag && $tag ne '';
    return '' unless $tag =~ /(?:FEMALE|MALE)/i;

    my $pair = uc $tag;
    $pair =~ s/_?AUTOSOME$//;
    $pair =~ s/_FEMALE/_SEX/;
    $pair =~ s/_MALE/_SEX/;
    return $pair;
}

sub beta {
    my ($row) = @_;
    if ($beta_col) {
        return numeric(value($row, $beta_col));
    }
    my $or = numeric(value($row, $or_col));
    return undef unless defined $or && $or > 0;
    return log($or);
}

sub value {
    my ($row, $col) = @_;
    return '' unless exists $idx{$col};
    return $row->[ $idx{$col} ] // '';
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
  perl diff_scz_female_male_gwas.pl [options]

Options:
  --input FILE.tsv.gz       Long-format sorted merged GWAS input
  --output FILE.tsv.gz      Female-vs-male differential output
  --manifest FILE.tsv       Run summary
  --beta-col NAME           Use an existing beta column. Default: log(OR)
  --or-col NAME             OR column used when --beta-col is absent. Default: OR
  --se-col NAME             Standard error column. Default: SE
  --rho FLOAT               Optional female/male effect correlation. Default: 0
  --limit N                 Write only the first N paired SNP results for testing

Formula:
  diff_beta = female_beta - male_beta
  diff_se   = sqrt(female_se^2 + male_se^2 - 2*rho*female_se*male_se)
  diff_z    = diff_beta / diff_se
  diff_p    = two-sided normal P value from diff_z
USAGE
}
