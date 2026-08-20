#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use JSON::PP qw(encode_json);
use Text::CSV;

my ($candidates, $population, $output, $summary);
my $input = '-';
my $snp_column = 'SNP';
my $min_r2 = 0.2;
my $candidate_proxies_only = 1;

GetOptions(
    'candidates=s'             => \$candidates,
    'population=s'             => \$population,
    'input=s'                  => \$input,
    'output=s'                 => \$output,
    'summary=s'                => \$summary,
    'snp-column=s'             => \$snp_column,
    'min-r2=f'                 => \$min_r2,
    'candidate-proxies-only!'  => \$candidate_proxies_only,
) or die usage();

die usage() unless defined $candidates && defined $population && defined $output;
$population = uc $population;
die "--population must be AFR, AMR, ASN, or EUR\n"
    unless $population =~ /\A(?:AFR|AMR|ASN|EUR)\z/;
die "--min-r2 must be at least 0.2 for the downloadable HaploReg v4 LD archives\n"
    unless $min_r2 >= 0.2 && $min_r2 <= 1;
$summary //= "$output.summary.json";

my %candidate = load_candidates($candidates, $snp_column);
die "No candidate rsIDs were read from $candidates\n" unless %candidate;

my $in_fh;
if ($input eq '-') {
    $in_fh = *STDIN;
}
elsif ($input =~ /\.gz\z/i) {
    require IO::Uncompress::Gunzip;
    $in_fh = IO::Uncompress::Gunzip->new($input)
        or die "Cannot gunzip $input\n";
}
else {
    open $in_fh, '<:raw', $input or die "Cannot read $input: $!\n";
}

open my $out_fh, '>:raw', $output or die "Cannot write $output: $!\n";
print {$out_fh} join("\t", qw(query_snp proxy_snp ld_population proxy_r2 proxy_dprime cache_min_r2)), "\n";

my ($archive_rows, $matched_queries, $proxy_edges, $candidate_edges) = (0, 0, 0, 0);
my %found;
while (my $line = <$in_fh>) {
    ++$archive_rows;
    $line =~ s/[\r\n]+\z//;
    my ($query, $proxies) = split /\t/, $line, 2;
    next unless defined $query;
    my $query_key = uc $query;
    next unless $candidate{$query_key};

    ++$matched_queries unless $found{$query_key}++;
    print {$out_fh} join("\t", $query_key, $query_key, $population, 1, 1, $min_r2), "\n";
    next unless defined $proxies && length $proxies;

    for my $edge (split /;/, $proxies) {
        my ($proxy, $r2, $dprime) = split /,/, $edge, 3;
        next unless defined $proxy && defined $r2 && $r2 =~ /\A(?:\d+(?:\.\d*)?|\.\d+)\z/;
        ++$proxy_edges;
        next if $r2 < $min_r2;
        my $proxy_key = uc $proxy;
        next if $candidate_proxies_only && !$candidate{$proxy_key};
        ++$candidate_edges;
        $dprime = '' unless defined $dprime;
        print {$out_fh} join("\t", $query_key, $proxy_key, $population, $r2, $dprime, $min_r2), "\n";
    }
}
close $out_fh or die "Cannot close $output: $!\n";
close $in_fh if $input ne '-';

my @missing = sort grep { !$found{$_} } keys %candidate;
my %stats = (
    population             => $population,
    archive_source         => $input,
    candidate_file         => $candidates,
    candidate_count        => scalar(keys %candidate),
    archive_rows_scanned   => $archive_rows,
    matched_query_count    => $matched_queries,
    missing_query_count    => scalar(@missing),
    missing_queries        => \@missing,
    archive_proxy_edges    => $proxy_edges,
    retained_candidate_edges => $candidate_edges,
    minimum_r2             => 0 + $min_r2,
    candidate_proxies_only => $candidate_proxies_only ? JSON::PP::true : JSON::PP::false,
    output                 => $output,
);
open my $summary_fh, '>:raw', $summary or die "Cannot write $summary: $!\n";
print {$summary_fh} JSON::PP->new->ascii->canonical->pretty->encode(\%stats);
close $summary_fh or die "Cannot close $summary: $!\n";

print "OUTPUT\t$output\n";
print "SUMMARY\t$summary\n";
print "CANDIDATES\t$stats{candidate_count}\n";
print "MATCHED_QUERIES\t$matched_queries\n";
print "MISSING_QUERIES\t", scalar(@missing), "\n";
print "RETAINED_CANDIDATE_EDGES\t$candidate_edges\n";

sub load_candidates {
    my ($path, $wanted_column) = @_;
    open my $fh, '<:raw', $path or die "Cannot read candidate file $path: $!\n";
    my $separator = $path =~ /\.csv\z/i ? ',' : "\t";
    my $csv = Text::CSV->new({ binary => 1, sep_char => $separator, auto_diag => 2 });
    my $header = $csv->getline($fh) or die "Candidate file is empty: $path\n";
    my %index = map { uc($header->[$_]) => $_ } 0 .. $#$header;
    my $column_key = uc $wanted_column;
    die "Candidate file $path lacks column $wanted_column\n" unless exists $index{$column_key};
    my %seen;
    while (my $row = $csv->getline($fh)) {
        my $value = $row->[$index{$column_key}];
        next unless defined $value;
        $value =~ s/^\s+|\s+$//g;
        next unless $value =~ /\Ars\d+\z/i;
        $seen{uc $value} = 1;
    }
    close $fh;
    return %seen;
}

sub usage {
    return <<'USAGE';
Usage:
  extract_haploreg_ld_for_candidates.pl \
    --candidates top_hit_candidates.csv --population EUR \
    --input LD_EUR.tsv.gz --output candidate_ld_eur.tsv [options]

The input may be '-' to read a decompressed HaploReg LD archive from STDIN.
The downloadable archives contain LD at r2>=0.2, so lower thresholds are
rejected. By default, only edges whose query and proxy are both candidates are
written; use --no-candidate-proxies-only to retain every proxy.
USAGE
}
