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
    '--title', 'Example & forest title',
    '--image-dir', $tmp,
    '--top-hits-csv', $csv,
);
is(system(@cmd), 0, 'forest HTML builder exits successfully');
ok(-s $html, 'self-contained forest HTML is nonempty');

open my $html_fh, '<:raw', $html or die "Cannot read $html: $!";
local $/;
my $content = <$html_fh>;
close $html_fh;

like($content, qr{<h1>Example &amp; forest title</h1>}, 'visible overall title is escaped');
like($content, qr{<figcaption>rs123 \(GENE1\)</figcaption>}, 'blank single-SNP label falls back to SNP and gene');
like($content, qr{src="data:image/png;base64,}, 'PNG is embedded as a data URI');
unlike($content, qr{src="forest_single_snp\.png"}, 'HTML does not depend on a relative sidecar path');

my ($embedded) = $content =~ m{data:image/png;base64,([^"]+)};
is(decode_base64($embedded // ''), $png_bytes, 'embedded PNG bytes match the source image');

done_testing();
