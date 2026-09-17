#!/usr/bin/env perl
use strict;
use warnings;

use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use File::Spec;
use Getopt::Long qw(GetOptions);

my ($queries, $pfile, $bfile, $plink2, $populations, $output);
my $min_r2 = 0;
my $window_kb = 1000;
my $quiet = 0;
GetOptions(
    'query-snps=s' => \$queries,
    'pfile=s' => \$pfile,
    'bfile=s' => \$bfile,
    'plink2=s' => \$plink2,
    'populations=s' => \$populations,
    'min-r2=f' => \$min_r2,
    'window-kb=f' => \$window_kb,
    'output=s' => \$output,
    'quiet!' => \$quiet,
) or die usage();

die usage() unless defined($queries) && length($queries) && defined($output)
    && (defined($pfile) xor defined($bfile));
my %seen;
my @queries;
for my $query (split /,/, $queries) {
    $query =~ s/^\s+|\s+$//g;
    next unless length($query) && !$seen{lc($query)}++;
    push @queries, $query;
}
die "--query-snps did not contain a query SNP\n" unless @queries;

my $resolver = File::Spec->catfile(dirname(__FILE__), 'resolve_plink2_local_ld.pl');
die "Missing PLINK2 LD resolver: $resolver\n" unless -f $resolver;
my $tmp = tempdir('plink2_multi_ld.XXXXXX', TMPDIR => 1, CLEANUP => 1);
my (@header, @rows, %header_idx, %proxy_owner);
for my $i (0 .. $#queries) {
    my $part = File::Spec->catfile($tmp, sprintf('query_%04d.tsv', $i + 1));
    my @cmd = (
        $^X, $resolver,
        '--query-snp', $queries[$i],
        '--min-r2', $min_r2,
        '--window-kb', $window_kb,
        '--output', $part,
        '--quiet',
    );
    push @cmd, ('--plink2', $plink2) if defined($plink2) && length($plink2);
    push @cmd, ('--populations', $populations) if defined($populations) && length($populations);
    push @cmd, defined($pfile) ? ('--pfile', $pfile) : ('--bfile', $bfile);
    system(@cmd) == 0 or die "PLINK2 LD calculation failed for $queries[$i]\n";
    open my $fh, '<', $part or die "Cannot read $part: $!\n";
    my $part_header = <$fh> // die "Empty LD cache for $queries[$i]\n";
    chomp $part_header;
    $part_header =~ s/\r$//;
    if (!@header) {
        @header = split /\t/, $part_header, -1;
        %header_idx = map { lc($header[$_]) => $_ } 0 .. $#header;
        die "LD cache is missing query_snp/proxy_snp columns\n"
            unless exists($header_idx{query_snp}) && exists($header_idx{proxy_snp});
    }
    die "Inconsistent LD cache header for $queries[$i]\n"
        unless $part_header eq join("\t", @header);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @fields = split /\t/, $line, -1;
        my $query = $fields[$header_idx{query_snp}] // '';
        my $proxy = $fields[$header_idx{proxy_snp}] // '';
        my $proxy_key = lc($proxy);
        if (length($proxy_key) && exists($proxy_owner{$proxy_key})
            && lc($proxy_owner{$proxy_key}) ne lc($query)) {
            die "Proxy $proxy occurs in more than one query window ($proxy_owner{$proxy_key}, $query). "
                . "Run the overlapping locus with one explicit LD reference SNP.\n";
        }
        $proxy_owner{$proxy_key} = $query if length $proxy_key;
        push @rows, $line;
    }
    close $fh;
}

make_path(dirname($output)) unless -d dirname($output);
my $temp_output = "$output.tmp.$$";
open my $out, '>', $temp_output or die "Cannot write $temp_output: $!\n";
print {$out} join("\t", @header), "\n", map { "$_\n" } @rows;
close $out or die "Cannot close $temp_output: $!\n";
rename $temp_output, $output or die "Cannot publish $output: $!\n";
print "LD_QUERY_COUNT\t", scalar(@queries), "\n" unless $quiet;
print "LD_CACHE_ROWS\t", scalar(@rows), "\n" unless $quiet;
print "LD_CACHE\t$output\n" unless $quiet;

sub usage {
    return <<'USAGE';
Usage: build_plink2_local_ld_cache.pl --query-snps rs1,rs2 (--pfile PREFIX | --bfile PREFIX) --output FILE [options]
  --plink2 EXE       PLINK 2 executable
  --populations LIST 1000 Genomes superpopulations, such as EUR
  --window-kb N      Maximum LD distance (default 1000)
  --min-r2 N         Minimum reported r2 (default 0)
  --quiet            Suppress summary output
USAGE
}
