#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use JSON::PP qw(encode_json);
use Test::More;

my $root = tempdir('oda_cleanup_XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $cleaner = "$Bin/cleanup_successful_oda_artifacts.pl";
my $png = "$root/plot.png";
my $html = "$root/plot.html";
write_file($png, 'PNG');
write_file($html, '<html>plot</html>');

my $stamp = '20261002_120000';
my $run = "$root/run_single_snp_with_gtf_$stamp";
my $upload = "$root/upload_single_snp_with_gtf_support_$stamp";
make_path($run, $upload);
write_file("$run/output.run.status.json", encode_json({state => 'completed', success => 1, complete => 1}));
write_file("$upload/upload.log", 'done');
is(run_cleaner($run, $upload), 0, 'completed run can clean matching helper folders');
ok(!-e $run && !-e $upload, 'run and upload folders are removed');
ok(-s $png && -s $html, 'final plot files are preserved');

$stamp = '20261002_120001';
$run = "$root/run_single_snp_with_gtf_$stamp";
$upload = "$root/upload_single_snp_with_gtf_support_$stamp";
make_path($run, $upload);
write_file("$run/output.run.status.json", encode_json({state => 'failed', success => 0, complete => 1}));
is(run_cleaner($run, $upload), 0, 'failed run does not turn cleanup into another error');
ok(-d $run && -d $upload, 'failed run diagnostics are preserved');

$stamp = '20261002_120002';
$run = "$root/run_manhattan_png_$stamp";
$upload = "$root/upload_manhattan_png_macro_$stamp";
make_path($run, $upload);
write_file("$run/output.run.status.json", encode_json({state => 'completed', success => 1, complete => 1}));
{
    local $ENV{KEEP_LOCAL_ODA_ARTIFACTS} = 1;
    is(run_cleaner($run, $upload), 0, 'explicit debug setting skips cleanup');
}
ok(-d $run && -d $upload, 'debug folders remain available');

my $nested = "$root/subdir/upload_manhattan_png_macro_$stamp";
make_path($nested);
isnt(run_cleaner($run, $nested), 0, 'nested path is rejected');
ok(-d $run && -d $nested, 'unsafe target does not remove the run or nested folder');
is(run_cleaner($run, $upload), 0, 'completed Manhattan run can be cleaned');
ok(!-e $run && !-e $upload, 'Manhattan run and upload folders are removed');

$stamp = '20261002_120003';
$run = "$root/run_local_hits_manhattan_png_$stamp";
$upload = "$root/upload_local_hits_support_$stamp";
my $autogen = "$root/.autogen_get_gtf_macro_local_mh_$stamp";
make_path($run, $upload, $autogen);
write_file("$run/output.run.status.json", encode_json({state => 'completed', success => 1, complete => 1}));
is(run_cleaner($run, $upload, $autogen), 0,
    'completed local Manhattan run cleans its generated GTF macro folder');
ok(!-e $run && !-e $upload && !-e $autogen,
    'local Manhattan run, upload, and generated macro folders are removed');

done_testing();

sub run_cleaner {
    my ($run_dir, @related) = @_;
    return system($^X, $cleaner,
        '--workdir', $root,
        '--run-dir', $run_dir,
        '--require-output', $png,
        '--require-output', $html,
        (map { ('--related-dir', $_) } @related));
}

sub write_file {
    my ($path, $contents) = @_;
    open my $fh, '>', $path or die "Cannot write $path: $!\n";
    print {$fh} $contents;
    close $fh or die "Cannot close $path: $!\n";
}
