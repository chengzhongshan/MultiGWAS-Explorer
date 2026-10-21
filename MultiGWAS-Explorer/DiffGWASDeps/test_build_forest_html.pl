#!/usr/bin/env perl
use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use MIME::Base64 qw(decode_base64);
use Test::More;

my $helper = File::Spec->catfile($Bin, 'build_forest_html.pl');
ok(-f $helper, 'forest HTML builder exists');

my $tmp = tempdir('forest_html_test_XXXX', TMPDIR => 1, CLEANUP => 1);
my $png = File::Spec->catfile($tmp, 'forest_single_snp.png');
my $manifest = File::Spec->catfile($tmp, 'forest.manifest.tsv');
my $csv = File::Spec->catfile($tmp, 'forest_top_hits.csv');
my $html = File::Spec->catfile($tmp, 'forest.html');

my $png_bytes = decode_base64(
    'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII='
);
open my $png_fh, '>:raw', $png or die "Cannot write $png: $!";
print {$png_fh} $png_bytes;
close $png_fh;

open my $manifest_fh, '>:encoding(UTF-8)', $manifest or die "Cannot write $manifest: $!";
print {$manifest_fh} "track_order\ttrack_id\ttrack_label\tpng_file\n";
print {$manifest_fh} " 1\tsingle_snp\t   \tforest_single_snp.png\n";
close $manifest_fh;

open my $csv_fh, '>:encoding(UTF-8)', $csv or die "Cannot write $csv: $!";
print {$csv_fh} "SNP,gene,comment\n";
print {$csv_fh} "rs123,GENE1,\"quoted, value\"\n";
close $csv_fh;

my @cmd = (
    $^X, $helper,
    '--manifest', $manifest,
    '--output', $html,
    '--title', 'Example GENE1 & GENE10 forest title',
    '--image-dir', $tmp,
    '--top-hits-csv', $csv,
);
is(system(@cmd), 0, 'forest HTML builder exits successfully');
ok(-s $html, 'self-contained forest HTML is nonempty');

open my $html_fh, '<:raw', $html or die "Cannot read $html: $!";
local $/;
my $content = <$html_fh>;
close $html_fh;

like($content, qr{<h1>Example <em class="gene-symbol">GENE1</em> &amp; GENE10 forest title</h1>},
    'known gene symbol is italicized in the visible overall title');
like($content, qr{<figcaption>rs123 \(<em class="gene-symbol">GENE1</em>\)</figcaption>},
    'blank single-SNP caption italicizes the fallback gene symbol');
like($content, qr{<title>Example GENE1 &amp; GENE10 forest title</title>},
    'document title metadata remains plain text');
unlike($content, qr{<em class="gene-symbol">GENE1</em>0},
    'gene matching does not italicize a substring of a longer symbol');
like($content, qr{alt="rs123 \(GENE1\) forest plot"},
    'image alternative text remains markup-free');
like($content, qr{src="data:image/png;base64,}, 'PNG is embedded as a data URI');
unlike($content, qr{src="forest_single_snp\.png"}, 'HTML does not depend on a relative sidecar path');

my ($embedded) = $content =~ m{data:image/png;base64,([^"]+)};
is(decode_base64($embedded // ''), $png_bytes, 'embedded PNG bytes match the source image');

my $forest_renderer = File::Spec->catfile($Bin, 'gnuplot', 'pdl_gnuplot_forest.pl');
my $render_csv = File::Spec->catfile($tmp, 'renderer.csv');
my $render_prefix = File::Spec->catfile($tmp, 'renderer');
open my $render_csv_fh, '>:encoding(UTF-8)', $render_csv or die "Cannot write $render_csv: $!";
print {$render_csv_fh} "SNP,gene,hit_class,BETA1,SE1,P1\n";
print {$render_csv_fh} "rs123,CYP26B1,DIFFERENTIAL,0.2,0.1,0.01\n";
close $render_csv_fh;

my $fake_gnuplot;
if ($^O =~ /MSWin32/i) {
    $fake_gnuplot = File::Spec->catfile($tmp, 'fake_gnuplot.cmd');
    open my $fake_fh, '>:raw', $fake_gnuplot or die "Cannot write $fake_gnuplot: $!";
    print {$fake_fh} "\@echo off\r\nexit /b 0\r\n";
    close $fake_fh;
}
else {
    $fake_gnuplot = File::Spec->catfile($tmp, 'fake_gnuplot.sh');
    open my $fake_fh, '>:raw', $fake_gnuplot or die "Cannot write $fake_gnuplot: $!";
    print {$fake_fh} "#!/bin/sh\nexit 0\n";
    close $fake_fh;
    chmod 0755, $fake_gnuplot or die "Cannot chmod $fake_gnuplot: $!";
}

my @render_cmd = (
    $^X, $forest_renderer,
    '--csv', $render_csv,
    '--out-prefix', $render_prefix,
    '--track-ids', 'TRACK1',
    '--track-labels', 'Track 1',
    '--track-beta-vars', 'BETA1',
    '--track-se-vars', 'SE1',
    '--track-p-vars', 'P1',
    '--gnuplot', $fake_gnuplot,
);
is(system(@render_cmd), 0, 'gnuplot forest renderer emits its script with a fake backend');
my $gp_content = slurp(File::Spec->catfile($tmp, 'renderer_single_snp.gp'));
like($gp_content, qr{set title "rs123 \(\{/:Italic CYP26B1\}\)" enhanced},
    'gnuplot forest title italicizes only the gene symbol');
my $renderer_manifest = slurp(File::Spec->catfile($tmp, 'renderer.manifest.tsv'));
like($renderer_manifest, qr{single_snp\trs123 \(CYP26B1\)\trenderer_single_snp\.png},
    'gnuplot manifest retains a plain-text SNP and gene label');

my $sas_macro = slurp(File::Spec->catfile($Bin, 'beta2OR_forest_plot.sas'));
like($sas_macro, qr{title font="Arial" "%superq\(plot_title_prefix\)" italic "%superq\(plot_title_gene\)" font="Arial" "%superq\(plot_title_suffix\)";},
    'SAS forest macro applies mixed normal and italic title segments');
my $sas_template = slurp(File::Spec->catfile($Bin, 'run_sas_oda_top_hits_forest_plot.sas'));
like($sas_template, qr{plot_title_gene=%superq\(forest_single_gene\)},
    'SAS forest template supplies the gene as a dedicated title segment');

done_testing();

sub slurp {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "Cannot read $path: $!";
    local $/;
    my $value = <$fh>;
    close $fh;
    return defined($value) ? $value : '';
}
