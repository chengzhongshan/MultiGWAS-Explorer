#!/usr/bin/env perl
use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IO::Compress::Gzip qw(gzip $GzipError);

my $gnuplot = $ENV{GNUPLOT} || 'gnuplot';
my $tmp = tempdir('multigwas_ld_heatmap_XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $plain = File::Spec->catfile($tmp, 'fixture.tsv');
my $gz = "$plain.gz";
my $gtf = File::Spec->catfile($tmp, 'fixture.gtf.tsv');
my $ld_cache = File::Spec->catfile($tmp, 'fixture.ld.tsv');
my $prefix = File::Spec->catfile($tmp, 'locus');

open my $fh, '>:raw', $plain or die "Cannot write $plain: $!\n";
print {$fh} "CHR\tBP\tSNP\tEUR_P\tEUR_Z\n";
print {$fh} "1\t100\trs100\t1e-6\t-4\n";
print {$fh} "1\t120\trs200\t1e-5\t2\n";
print {$fh} "1\t140\trs300\t1e-4\t4\n";
close $fh or die "Cannot close $plain: $!\n";
gzip $plain => $gz or die "Cannot gzip fixture: $GzipError\n";

open my $gf, '>:raw', $gtf or die "Cannot write $gtf: $!\n";
print {$gf} "chr\tgenesymbol\tgene\tst\ten\ttype\n";
print {$gf} "1\tGENE1\tGENE1\t80\t160\tgene\n";
print {$gf} "1\tGENE1\tGENE1\t90\t110\texon\n";
close $gf or die "Cannot close $gtf: $!\n";

open my $lf, '>:raw', $ld_cache or die "Cannot write $ld_cache: $!\n";
print {$lf} "query_snp\tproxy_snp\tld_population\tproxy_r2\n";
print {$lf} "RS100\tRS200\tEUR\t0.55\n";
print {$lf} "RS100\tRS300\tEUR\t0.91\n";
close $lf or die "Cannot close $ld_cache: $!\n";
my $resolver = File::Spec->catfile($Bin, 'resolve_haploreg_high_ld.pl');
open my $resolver_fh, '-|', $^X, $resolver,
    '--query-snps', 'rs100', '--population', 'EUR', '--min-r2', '0.5',
    '--local-cache', $ld_cache, '--no-web-fallback'
    or die "Cannot start LD resolver: $!\n";
my $resolver_output = do { local $/; <$resolver_fh> };
close $resolver_fh or die "LD resolver fixture failed\n";
die "LD resolver did not preserve numeric r2 pairs\n"
    unless $resolver_output =~ /^LD_R2_PAIRS\tRS300:0\.91,RS200:0\.55$/m;

my $renderer = File::Spec->catfile($Bin, 'gnuplot', 'pdl_gunplot_local_locus.pl');
my @cmd = (
    $^X, $renderer,
    '--data', $gz,
    '--snp', 'rs100',
    '--label-snps', 'rs100',
    '--out-prefix', $prefix,
    '--window-bp', 100,
    '--pcols', 'EUR_P',
    '--zcols', 'EUR_Z',
    '--labels', 'EUR',
    '--gtf', $gtf,
    '--gnuplot', $gnuplot,
    '--ld-snps', 'rs200,rs300',
    '--ld-r2-values', 'rs200:0.55,rs300:0.91',
    '--ld-display-mode', 'heatmap',
    '--ld-population', 'EUR',
);
system(@cmd) == 0 or die "LD heatmap renderer failed (exit " . ($? >> 8) . ")\n";

for my $path ("$prefix.png", "$prefix.gp", "$prefix.plot.tsv", "$prefix.manifest.tsv") {
    die "Expected output is missing or empty: $path\n" unless -s $path;
}
open my $gp, '<:raw', "$prefix.gp" or die $!;
my $gp_text = do { local $/; <$gp> };
close $gp;
die "Association Z-score colorbar is missing\n" unless $gp_text =~ /set cblabel 'Z score'/;
die "Separate LD inset title is missing\n" unless $gp_text =~ /LD r\^2 \(EUR\)/;
die "LD RGB-variable overlay is missing\n" unless $gp_text =~ /using \(\(\$9>=0\)\?\$1:1\/0\):2:10/;
my ($ld_low_color) = $gp_text =~ /set object 8000 .*?fc rgb '(#[0-9a-f]{6})'/;
die "LD inset palette was not emitted\n" unless defined $ld_low_color;
die "LD palette unexpectedly reused an association endpoint color\n"
    if $ld_low_color =~ /^(?:#63d67f|#63d8d2|#ffbf00|#ff5b00|#df1f2d)$/;

open my $pt, '<:raw', "$prefix.plot.tsv" or die $!;
my $header = <$pt> // '';
close $pt;
die "LD numeric/color columns are missing\n" unless $header =~ /\tLD_R2\tLD_RGB\s*$/;

print "Optional LD heatmap rendering: PASS\n";
