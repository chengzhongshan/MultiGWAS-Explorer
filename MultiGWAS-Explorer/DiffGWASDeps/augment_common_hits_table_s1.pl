#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Text::ParseWords qw(parse_line);

sub usage {
    return <<"USAGE";
Usage:
  perl DiffGWASDeps/augment_common_hits_table_s1.pl [options]

Options:
  --input-csv FILE     Existing common-hit CSV. Default:
                       manuscript_assets/tables/Table_S1_all_common_association_loci.csv
  --wide FILE          Wide standardized GWAS TSV.GZ. Required.
  --output-csv FILE    Output augmented CSV. Required.
USAGE
}

my %opt = (
    input_csv => 'manuscript_assets/tables/Table_S1_all_common_association_loci.csv',
);

GetOptions(
    'input-csv=s'  => \$opt{input_csv},
    'wide=s'       => \$opt{wide},
    'output-csv=s' => \$opt{output_csv},
) or die usage();

die usage() unless $opt{wide} && $opt{output_csv};
die "Input CSV not found: $opt{input_csv}\n" unless -s $opt{input_csv};
die "Wide file not found: $opt{wide}\n" unless -s $opt{wide};

my ($header, $rows) = read_csv_rows($opt{input_csv});
my %by_snp = map { uc(trim($_->{SNP} // '')) => $_ } @$rows;
die "No SNP rows loaded from $opt{input_csv}\n" unless %by_snp;

my @wanted_cols = qw(
  ALL_GROUP1_BETA ALL_GROUP1_SE ALL_GROUP1_P
  ALL_GROUP2_BETA ALL_GROUP2_SE ALL_GROUP2_P
  ALL_DIFF_BETA ALL_DIFF_SE ALL_DIFF_P ALL_STD_DIFF_P
  EUR_GROUP1_BETA EUR_GROUP1_SE EUR_GROUP1_P
  EUR_GROUP2_BETA EUR_GROUP2_SE EUR_GROUP2_P
  EUR_DIFF_BETA EUR_DIFF_SE EUR_DIFF_P EUR_STD_DIFF_P
  ASN_GROUP1_BETA ASN_GROUP1_SE ASN_GROUP1_P
  ASN_GROUP2_BETA ASN_GROUP2_SE ASN_GROUP2_P
  ASN_DIFF_BETA ASN_DIFF_SE ASN_DIFF_P ASN_STD_DIFF_P
);

my $z = IO::Uncompress::Gunzip->new($opt{wide})
  or die "Cannot open $opt{wide}: $GunzipError\n";
my $wide_header = <$z>;
die "Wide file is empty: $opt{wide}\n" unless defined $wide_header;
chomp $wide_header;
$wide_header =~ s/\r$//;
my @wide_cols = split /\t/, $wide_header, -1;
my %widx = map { $wide_cols[$_] => $_ } 0 .. $#wide_cols;
for my $c (qw(SNP CHR BP), @wanted_cols) {
    die "Wide file is missing required column $c\n" unless exists $widx{$c};
}

my %wide_by_snp;
while (my $line = <$z>) {
    chomp $line;
    $line =~ s/\r$//;
    next unless length $line;
    my @f = split /\t/, $line, -1;
    my $snp = uc(trim($f[$widx{SNP}] // ''));
    next unless length $snp && exists $by_snp{$snp};
    my %row;
    for my $c (@wanted_cols, qw(CHR BP SNP)) {
        $row{$c} = $f[$widx{$c}];
    }
    $wide_by_snp{$snp} = \%row;
}
close $z;

my @out_header = (
    qw(CHR BP SNP gene Smallest_ASSOC_P),
    qw(ALL_STD_DIFF_P EUR_STD_DIFF_P ASN_STD_DIFF_P),
    qw(ALL_FEMALE_BETA ALL_FEMALE_SE ALL_FEMALE_P),
    qw(ALL_MALE_BETA ALL_MALE_SE ALL_MALE_P),
    qw(ALL_DIFF_BETA ALL_DIFF_SE ALL_DIFF_P),
    qw(EUR_FEMALE_BETA EUR_FEMALE_SE EUR_FEMALE_P),
    qw(EUR_MALE_BETA EUR_MALE_SE EUR_MALE_P),
    qw(EUR_DIFF_BETA EUR_DIFF_SE EUR_DIFF_P),
    qw(ASN_FEMALE_BETA ASN_FEMALE_SE ASN_FEMALE_P),
    qw(ASN_MALE_BETA ASN_MALE_SE ASN_MALE_P),
    qw(ASN_DIFF_BETA ASN_DIFF_SE ASN_DIFF_P),
);

open my $out, '>', $opt{output_csv} or die "Cannot write $opt{output_csv}: $!\n";
print {$out} join(",", map { csv_escape($_) } @out_header), "\n";
for my $src (@$rows) {
    my $snp = uc(trim($src->{SNP} // ''));
    my $wide = $wide_by_snp{$snp}
      or die "Could not find SNP $snp in wide file while building Table S1\n";
    my %out_row = (
        CHR              => value_or($src->{CHR}, $wide->{CHR}),
        BP               => value_or($src->{BP}, $wide->{BP}),
        SNP              => value_or($src->{SNP}, $wide->{SNP}),
        gene             => value_or($src->{gene}, ''),
        Smallest_ASSOC_P => value_or($src->{COMMON_ASSOC_P}, value_or($src->{focus_signal}, '')),
        ALL_STD_DIFF_P   => value_or($wide->{ALL_STD_DIFF_P}, value_or($src->{ALL_STD_DIFF_P}, '')),
        EUR_STD_DIFF_P   => value_or($wide->{EUR_STD_DIFF_P}, ''),
        ASN_STD_DIFF_P   => value_or($wide->{ASN_STD_DIFF_P}, ''),
        ALL_FEMALE_BETA  => value_or($wide->{ALL_GROUP1_BETA}, ''),
        ALL_FEMALE_SE    => value_or($wide->{ALL_GROUP1_SE}, ''),
        ALL_FEMALE_P     => value_or($wide->{ALL_GROUP1_P}, ''),
        ALL_MALE_BETA    => value_or($wide->{ALL_GROUP2_BETA}, ''),
        ALL_MALE_SE      => value_or($wide->{ALL_GROUP2_SE}, ''),
        ALL_MALE_P       => value_or($wide->{ALL_GROUP2_P}, ''),
        ALL_DIFF_BETA    => value_or($wide->{ALL_DIFF_BETA}, ''),
        ALL_DIFF_SE      => value_or($wide->{ALL_DIFF_SE}, ''),
        ALL_DIFF_P       => value_or($wide->{ALL_DIFF_P}, ''),
        EUR_FEMALE_BETA  => value_or($wide->{EUR_GROUP1_BETA}, ''),
        EUR_FEMALE_SE    => value_or($wide->{EUR_GROUP1_SE}, ''),
        EUR_FEMALE_P     => value_or($wide->{EUR_GROUP1_P}, ''),
        EUR_MALE_BETA    => value_or($wide->{EUR_GROUP2_BETA}, ''),
        EUR_MALE_SE      => value_or($wide->{EUR_GROUP2_SE}, ''),
        EUR_MALE_P       => value_or($wide->{EUR_GROUP2_P}, ''),
        EUR_DIFF_BETA    => value_or($wide->{EUR_DIFF_BETA}, ''),
        EUR_DIFF_SE      => value_or($wide->{EUR_DIFF_SE}, ''),
        EUR_DIFF_P       => value_or($wide->{EUR_DIFF_P}, ''),
        ASN_FEMALE_BETA  => value_or($wide->{ASN_GROUP1_BETA}, ''),
        ASN_FEMALE_SE    => value_or($wide->{ASN_GROUP1_SE}, ''),
        ASN_FEMALE_P     => value_or($wide->{ASN_GROUP1_P}, ''),
        ASN_MALE_BETA    => value_or($wide->{ASN_GROUP2_BETA}, ''),
        ASN_MALE_SE      => value_or($wide->{ASN_GROUP2_SE}, ''),
        ASN_MALE_P       => value_or($wide->{ASN_GROUP2_P}, ''),
        ASN_DIFF_BETA    => value_or($wide->{ASN_DIFF_BETA}, ''),
        ASN_DIFF_SE      => value_or($wide->{ASN_DIFF_SE}, ''),
        ASN_DIFF_P       => value_or($wide->{ASN_DIFF_P}, ''),
    );
    print {$out} join(",", map { csv_escape($out_row{$_}) } @out_header), "\n";
}
close $out;

print "Wrote augmented Table S1 CSV: $opt{output_csv}\n";
print "Rows: " . scalar(@$rows) . "\n";

sub read_csv_rows {
    my ($file) = @_;
    open my $fh, '<', $file or die "Cannot read $file: $!\n";
    my $header_line = <$fh>;
    die "CSV is empty: $file\n" unless defined $header_line;
    chomp $header_line;
    $header_line =~ s/\r$//;
    my @header = map { trim($_) } parse_line(',', 0, $header_line);
    my @rows;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @vals = parse_line(',', 0, $line);
        my %row;
        @row{@header} = @vals;
        push @rows, \%row;
    }
    close $fh;
    return (\@header, \@rows);
}

sub csv_escape {
    my ($v) = @_;
    $v = '' unless defined $v;
    if ($v =~ /[",\n]/) {
        $v =~ s/"/""/g;
        return qq{"$v"};
    }
    return $v;
}

sub trim {
    my ($x) = @_;
    $x = '' unless defined $x;
    $x =~ s/^\s+//;
    $x =~ s/\s+$//;
    return $x;
}

sub value_or {
    my ($a, $b) = @_;
    return $a if defined $a && $a ne '';
    return $b if defined $b;
    return '';
}
