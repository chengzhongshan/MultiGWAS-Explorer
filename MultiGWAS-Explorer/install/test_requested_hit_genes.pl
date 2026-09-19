#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use Text::CSV;
use Test::More;

my $dir = tempdir(CLEANUP => 1);
my $wide = "CHR\tBP\tSNP\tEUR_DIFF_P\n"
    . "1\t150\trs1\t0.5\n1\t160\trs2\t0.5\n"
    . "1\t350\trs3\t0.5\n23\t150\trs4\t0.5\n2\t150\trs5\t0.5\n";
my $gtf = join('', map {
    my ($chr, $start, $end, $gene) = @$_;
    join("\t", $chr, 'fixture', 'gene', $start, $end, '.', '+', '.',
        qq{gene_id "$gene"; gene_name "$gene"; gene_type "protein_coding";}) . "\n"
} ['chr1',100,200,'GENE1'], ['chr1',400,500,'GENE2'], ['chrX',100,200,'XGENE']);
gzip \$wide => "$dir/wide.tsv.gz" or die $GzipError;
gzip \$gtf => "$dir/genes.gtf.gz" or die $GzipError;

sub export_hits {
    my ($name, @extra) = @_;
    my $path = "$dir/$name.csv";
    my @cmd = ($^X, "$Bin/../DiffGWASDeps/generate_requested_top_hits_csv.pl",
        '--input', "$dir/wide.tsv.gz", '--output', $path,
        '--gene-annotation-gtf', "$dir/genes.gtf.gz",
        '--target-snps', 'rs1,rs2,rs3,rs4,rs5', '--top-hit-focus-pvar', 'EUR_DIFF_P', @extra);
    is(system(@cmd), 0, "$name export succeeds") or return {};
    open my $fh, '<', $path or die $!;
    my $csv = Text::CSV->new({binary => 1});
    $csv->column_names($csv->getline($fh));
    my %rows;
    while (my $row = $csv->getline_hr($fh)) { $rows{$row->{SNP}} = $row; }
    close $fh;
    return \%rows;
}

my $rows = export_hits('annotated');
is(scalar keys %$rows, 5, 'explicit targets are retained regardless of significance');
is($rows->{rs1}{gene}, 'GENE1', 'overlapping gene is annotated');
is($rows->{rs2}{gene}, 'GENE1', 'multiple SNPs retain the same gene');
is($rows->{rs3}{gene}, 'GENE2', 'intergenic SNP gets nearest protein-coding gene');
is($rows->{rs4}{gene}, 'XGENE', 'chromosome 23 matches chrX annotation');
is($rows->{rs5}{gene}, '', 'absent chromosome does not borrow another chromosome gene');
is($rows->{rs1}{gene_source}, 'GTF', 'GTF provenance is preserved');
is($rows->{rs3}{snp_gene}, 'rs3:GENE2', 'combined label uses the resolved gene');
my $overridden = export_hits('overridden', '--target-snp-genes', 'rs1:CUSTOM');
is($overridden->{rs1}{gene}, 'CUSTOM', 'explicit gene override takes priority');
is($overridden->{rs2}{gene}, 'GENE1', 'other requested genes are still annotated');
done_testing;
