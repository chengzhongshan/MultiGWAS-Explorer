#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use FindBin qw($Bin);
use Cwd qw(abs_path);

my $cleaner = abs_path("$Bin/../clean_pipeline_temporary_files.pl");
my $tmp = tempdir(CLEANUP => 1);
chdir $tmp or die $!;
sub put {
    my ($path) = @_;
    open my $fh, '>', $path or die "$path: $!";
    print {$fh} "fixture\n";
    close $fh;
}
sub run_cleaner {
    open my $fh, '-|', $^X, $cleaner, @_ or die $!;
    my $out = do { local $/; <$fh> };
    close $fh;
    is($? >> 8, 0, 'cleanup exits successfully');
    return $out;
}
my $old = 'run_sas_oda_local_top_hits_with_gtf.20260901_120000.part1.sas';
my $fresh = 'sas_submit_result_ab12.json';
my $tracked = 'sas_inline_code_ab12.sas';
my $stage = 'upload_top_hits_batches_20260901_120000';
my $unknown = 'upload_manhattan_subset_20260901_120000';
make_path($stage, $unknown, 'cache', 'run_local_hits_with_gtf_20260901_120000');
put($_) for ($old, $fresh, $tracked, "$stage/input.csv", "$unknown/result.png",
    'final.png', 'final.html', 'cache/ld.tsv', 'run_local_hits_with_gtf_20260901_120000/output.html');
utime(time - 172800, time - 172800, $old);
system('git', 'init', '-q') == 0 or die 'git init failed';
system('git', 'add', $tracked) == 0 or die 'git add failed';
my $preview = run_cleaner();
like($preview, qr/WOULD DELETE\t\Q$old\E/, 'old generated script is previewed');
ok(-f $old, 'preview preserves files');
unlike($preview, qr/WOULD DELETE\t\Q$fresh\E/, 'recent file is protected by default');
run_cleaner('--apply');
ok(!-e $old && -f $fresh, 'default application removes only old temporary files');
my $link = 'target_snp_augmented_local_gtf_20260901_120000.tsv';
my $linked = symlink('final.png', $link);
run_cleaner('--min-age-hours', '0', '--apply');
ok(!-e $fresh && !-e $stage, 'explicit zero age removes recent helpers and staging');
ok(-f $tracked, 'tracked generated-looking source is preserved');
ok(-f "$unknown/result.png", 'unknown staging contents are preserved');
ok(-f $_, "preserved $_") for ('final.png', 'final.html', 'cache/ld.tsv',
    'run_local_hits_with_gtf_20260901_120000/output.html');
SKIP: {
    skip 'symlinks unavailable', 1 unless $linked;
    ok(-l $link && -f 'final.png', 'symlink and target preserved');
}
done_testing();
