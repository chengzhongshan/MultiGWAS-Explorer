#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);

sub usage {
    return <<'USAGE';
Usage:
  perl DiffGWASDeps/convert_ld_pruned_leads_to_common_loci_csv.pl \
    --input-tsv benchmark/reviewer_pgc_common_ld_cached_r2ge0p2_leads.tsv \
    --output-csv /path/to/Supplementary_Table_S12.csv

Converts a validated LD-pruned common-association lead TSV into the formatted
Table S12 column schema. The output retains the original association blocks and
adds LD-selection and MAF provenance fields.
USAGE
}

my %opt;
GetOptions(
    'input-tsv=s'  => \$opt{input_tsv},
    'output-csv=s' => \$opt{output_csv},
) or die usage();

die usage() unless $opt{input_tsv} && $opt{output_csv};
die "Input TSV not found: $opt{input_tsv}\n" unless -s $opt{input_tsv};

my @metadata_cols = qw(
  LD_LEAD_RANK INDEPENDENCE_METHOD LD_POPULATIONS LD_POPULATION_RULE
  LD_QUERY_STATUS LD_R2_THRESHOLD LD_FALLBACK_BP selected_maf maf_source
);
my @association_cols = qw(
  CHR BP SNP gene Smallest_ASSOC_P
  ALL_STD_DIFF_P EUR_STD_DIFF_P ASN_STD_DIFF_P
  ALL_FEMALE_BETA ALL_FEMALE_SE ALL_FEMALE_P
  ALL_MALE_BETA ALL_MALE_SE ALL_MALE_P
  ALL_DIFF_BETA ALL_DIFF_SE ALL_DIFF_P
  EUR_FEMALE_BETA EUR_FEMALE_SE EUR_FEMALE_P
  EUR_MALE_BETA EUR_MALE_SE EUR_MALE_P
  EUR_DIFF_BETA EUR_DIFF_SE EUR_DIFF_P
  ASN_FEMALE_BETA ASN_FEMALE_SE ASN_FEMALE_P
  ASN_MALE_BETA ASN_MALE_SE ASN_MALE_P
  ASN_DIFF_BETA ASN_DIFF_SE ASN_DIFF_P
);
my @output_cols = (@metadata_cols, @association_cols);

my %source_for = (
    Smallest_ASSOC_P => 'COMMON_ASSOC_P',
    ALL_FEMALE_BETA  => 'ALL_GROUP1_BETA',
    ALL_FEMALE_SE    => 'ALL_GROUP1_SE',
    ALL_FEMALE_P     => 'ALL_GROUP1_P',
    ALL_MALE_BETA    => 'ALL_GROUP2_BETA',
    ALL_MALE_SE      => 'ALL_GROUP2_SE',
    ALL_MALE_P       => 'ALL_GROUP2_P',
    EUR_FEMALE_BETA  => 'EUR_GROUP1_BETA',
    EUR_FEMALE_SE    => 'EUR_GROUP1_SE',
    EUR_FEMALE_P     => 'EUR_GROUP1_P',
    EUR_MALE_BETA    => 'EUR_GROUP2_BETA',
    EUR_MALE_SE      => 'EUR_GROUP2_SE',
    EUR_MALE_P       => 'EUR_GROUP2_P',
    ASN_FEMALE_BETA  => 'ASN_GROUP1_BETA',
    ASN_FEMALE_SE    => 'ASN_GROUP1_SE',
    ASN_FEMALE_P     => 'ASN_GROUP1_P',
    ASN_MALE_BETA    => 'ASN_GROUP2_BETA',
    ASN_MALE_SE      => 'ASN_GROUP2_SE',
    ASN_MALE_P       => 'ASN_GROUP2_P',
);

open my $in, '<', $opt{input_tsv}
  or die "Cannot read $opt{input_tsv}: $!\n";
my $header_line = <$in>;
die "Input TSV is empty: $opt{input_tsv}\n" unless defined $header_line;
chomp $header_line;
$header_line =~ s/\r$//;
my @header = split /\t/, $header_line, -1;
my %index = map { $header[$_] => $_ } 0 .. $#header;

my %required = map { $_ => 1 } @metadata_cols;
$required{maf_filter_decision} = 1;
for my $out_col (@association_cols) {
    my $source = $source_for{$out_col} // $out_col;
    $required{$source} = 1;
}
for my $name (sort keys %required) {
    die "Input TSV is missing required column $name\n"
      unless exists $index{$name};
}

my @rows;
my (%seen_snp, %seen_rank);
while (my $line = <$in>) {
    chomp $line;
    $line =~ s/\r$//;
    next unless length $line;
    my @field = split /\t/, $line, -1;
    my %source;
    @source{@header} = @field;

    my $snp = trim($source{SNP});
    my $rank = trim($source{LD_LEAD_RANK});
    die "Encountered a lead row without SNP\n" unless length $snp;
    die "Encountered non-numeric LD_LEAD_RANK for $snp: $rank\n"
      unless $rank =~ /^\d+$/ && $rank > 0;
    die "Duplicate SNP in lead input: $snp\n" if $seen_snp{uc $snp}++;
    die "Duplicate LD lead rank in input: $rank\n" if $seen_rank{$rank}++;
    die "LD lead $snp did not pass the MAF safeguard\n"
      unless uc(trim($source{maf_filter_decision})) eq 'PASS';

    my %out;
    for my $out_col (@output_cols) {
        my $source_col = $source_for{$out_col} // $out_col;
        $out{$out_col} = defined $source{$source_col} ? $source{$source_col} : '';
    }
    push @rows, \%out;
}
close $in;

die "No LD-pruned lead rows were loaded\n" unless @rows;
@rows = sort { $a->{LD_LEAD_RANK} <=> $b->{LD_LEAD_RANK} } @rows;
for my $i (0 .. $#rows) {
    my $expected = $i + 1;
    die "LD ranks are not continuous at expected rank $expected\n"
      unless $rows[$i]{LD_LEAD_RANK} == $expected;
}

open my $out, '>', $opt{output_csv}
  or die "Cannot write $opt{output_csv}: $!\n";
print {$out} join(',', map { csv_escape($_) } @output_cols), "\n";
for my $row (@rows) {
    print {$out} join(',', map { csv_escape($row->{$_}) } @output_cols), "\n";
}
close $out;

print "Wrote LD-pruned common-association CSV: $opt{output_csv}\n";
print "Rows: " . scalar(@rows) . "\n";
print "LD rank range: 1-" . scalar(@rows) . "\n";

sub csv_escape {
    my ($value) = @_;
    $value = '' unless defined $value;
    if ($value =~ /[",\r\n]/) {
        $value =~ s/"/""/g;
        return qq{"$value"};
    }
    return $value;
}

sub trim {
    my ($value) = @_;
    $value = '' unless defined $value;
    $value =~ s/^\s+//;
    $value =~ s/\s+$//;
    return $value;
}
