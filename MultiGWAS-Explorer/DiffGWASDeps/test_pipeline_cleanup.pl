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
my $action_old = 'sas_action_result_cd34.json';
my $action_fresh = 'sas_action_runner_ef56.py';
my $tracked = 'sas_inline_code_ab12.sas';
my $stage = 'upload_top_hits_batches_20260901_120000';
my $unknown = 'upload_manhattan_subset_20260901_120000';
my $old_run = 'run_local_hits_with_gtf_20260901_120000';
my $old_part = 'run_local_hits_with_gtf_20260901_120000_part1';
my $fresh_run = 'run_local_hits_with_gtf_20260902_120000';
my $tracked_run = 'run_local_hits_with_gtf_20260903_120000';
make_path($stage, $unknown, 'cache', "$old_run/nested", $old_part, $fresh_run,
    $tracked_run, 'run_local_hits_with_gtf_not_a_timestamp');
put($_) for ($old, $fresh, $action_old, $action_fresh, $tracked,
    "$stage/input.csv", "$unknown/result.png",
    'final.png', 'final.html', 'cache/ld.tsv', "$old_run/output.html.info.txt",
    "$old_run/nested/output.run.status.json", "$old_part/output.html.info.txt",
    "$fresh_run/output.html.info.txt", "$tracked_run/output.run.status.json",
    'run_local_hits_with_gtf_not_a_timestamp/keep.txt');
my $old_time = time - 172800;
utime($old_time, $old_time, $old, $action_old, "$old_run/output.html.info.txt",
    "$old_run/nested/output.run.status.json", "$old_part/output.html.info.txt",
    "$tracked_run/output.run.status.json");
utime($old_time, $old_time, "$old_run/nested", $old_run, $old_part, $tracked_run);
system('git', 'init', '-q') == 0 or die 'git init failed';
system('git', 'add', $tracked, "$tracked_run/output.run.status.json") == 0
    or die 'git add failed';
my $preview = run_cleaner();
like($preview, qr/WOULD DELETE\t\Q$old\E/, 'old generated script is previewed');
like($preview, qr/WOULD DELETE\t\Q$action_old\E/, 'old SAS action result is previewed');
like($preview, qr/WOULD RMDIR\t\Q$old_run\E/, 'old local-GTF run directory is previewed');
ok(-f $old, 'preview preserves files');
unlike($preview, qr/WOULD DELETE\t\Q$fresh\E/, 'recent file is protected by default');
run_cleaner('--apply');
ok(!-e $old && !-e $action_old && -f $fresh && -f $action_fresh,
    'default application removes only old temporary files');
ok(!-e $old_run && !-e $old_part, 'old local-GTF run directories are removed recursively');
ok(-d $fresh_run, 'recent local-GTF run directory is protected by default');
my $link = 'target_snp_augmented_local_gtf_20260901_120000.tsv';
my $linked = symlink('final.png', $link);
run_cleaner('--min-age-hours', '0', '--apply');
ok(!-e $fresh && !-e $action_fresh && !-e $stage && !-e $fresh_run,
    'explicit zero age removes recent helpers, SAS actions, staging and run directories');
ok(-f $tracked, 'tracked generated-looking source is preserved');
ok(-f "$tracked_run/output.run.status.json", 'run directory containing a tracked file is preserved');
ok(-f "$unknown/result.png", 'unknown staging contents are preserved');
ok(-f $_, "preserved $_") for ('final.png', 'final.html', 'cache/ld.tsv',
    'run_local_hits_with_gtf_not_a_timestamp/keep.txt');
SKIP: {
    skip 'symlinks unavailable', 1 unless $linked;
    ok(-l $link && -f 'final.png', 'symlink and target preserved');
}
done_testing();
