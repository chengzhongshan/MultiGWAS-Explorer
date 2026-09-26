#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename qw(dirname);
use File::Spec;
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use Test::More;

my $ancestor = $Bin;
my $have_tools = 0;
for (1 .. 4) {
    $ancestor = dirname($ancestor);
    my $bin = File::Spec->catdir($ancestor, 'local', 'bin');
    if ((-x File::Spec->catfile($bin, 'tabix') || -x File::Spec->catfile($bin, 'tabix.exe'))
        && (-x File::Spec->catfile($bin, 'bgzip') || -x File::Spec->catfile($bin, 'bgzip.exe'))) {
        $have_tools = 1;
        last;
    }
}
$have_tools = 1 if -x '/usr/bin/tabix' && -x '/usr/bin/bgzip';
plan skip_all => 'native tabix and bgzip are not installed' unless $have_tools;

my $dir = tempdir('gtf tools XXXXX', TMPDIR => 1, CLEANUP => 1);
my $source = "$dir/tiny.gtf.gz";
my $output = "$dir/region.tsv";
my $gtf = join("\t", 'chr1', 'test', 'gene', 100, 200, '.', '+', '.',
    'gene_id "GENE1"; gene_name "TEST1"; gene_type "protein_coding";') . "\n";
gzip \$gtf => $source or die $GzipError;
is(system($^X, "$Bin/../DiffGWASDeps/extract_gencode_gtf_subset.pl",
    '--gtf-gz', $source, '--output', $output,
    '--region', '1:110:190', '--reference-build', 'hg38'), 0,
    'GTF extraction discovers native tabix and bgzip without explicit paths');
ok(-s $output, 'GTF region output is nonempty');
my $indexed = "$dir/tiny.gtf.bgz";
ok(-s $indexed && -s "$indexed.tbi", 'native tools build a reusable BGZF GTF index');
done_testing();
