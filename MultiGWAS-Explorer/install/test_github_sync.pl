#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Copy qw(copy);
use Cwd qw(getcwd abs_path);
use IPC::Open3;
use IO::Select;
use Symbol qw(gensym);
use Test::More;

# Exercise real Git synchronization using two disposable clones and a bare
# local remote. No GitHub connection or real user repository is modified.
my $dir = tempdir('multigwas_git_sync_XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $remote = "$dir/remote.git";
my $local = "$dir/local copy";
my $other = "$dir/other computer";
my $script = abs_path("$Bin/../../update_github_upload.sh");
local $ENV{GIT_CONFIG_GLOBAL} = "$dir/empty.gitconfig";
local $ENV{GIT_CONFIG_NOSYSTEM} = 1;
local $ENV{GIT_TERMINAL_PROMPT} = 0;
local $ENV{REMOTE_NAME} = 'origin';
local $ENV{BRANCH_NAME} = 'main';
delete local $ENV{GIT_INDEX_FILE};

sub run_in {
    my ($cwd, @cmd) = @_;
    my $previous = getcwd();
    chdir $cwd or die "chdir $cwd: $!";
    my $err = gensym;
    my $pid = open3(my $in, my $out, $err, @cmd);
    close $in;
    chdir $previous or die "chdir $previous: $!";
    my $readers = IO::Select->new($out, $err);
    my $output = '';
    while ($readers->count) {
        for my $fh ($readers->can_read) {
            my $n = sysread($fh, my $buffer, 8192);
            die "read child output: $!" unless defined $n;
            if ($n) { $output .= $buffer; }
            else { $readers->remove($fh); close $fh; }
        }
    }
    waitpid($pid, 0);
    return ($? >> 8, $output);
}

sub command {
    my ($cwd, @cmd) = @_;
    my ($status, $output) = run_in($cwd, @cmd);
    die "@cmd failed ($status):\n$output" if $status;
    $output =~ s/\s+$//;
    return $output;
}

sub write_file {
    my ($path, $text) = @_;
    open my $fh, '>:raw', $path or die "write $path: $!";
    print {$fh} $text;
    close $fh or die "close $path: $!";
}

sub read_file {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "read $path: $!";
    return do { local $/; <$fh> };
}

sub commit_all {
    my ($cwd, $message) = @_;
    command($cwd, 'git', 'add', '-A');
    command($cwd, 'git', 'commit', '-qm', $message);
}

sub remote_edit {
    my ($path, $text) = @_;
    write_file("$other/$path", $text);
    commit_all($other, "Other computer changes $path");
    command($other, 'git', 'push', '-q', 'origin', 'main');
}

sub synchronize {
    my (@args) = @_;
    return run_in($dir, '/bin/bash', "$local/update_github_upload.sh", @args);
}

command($dir, 'git', 'init', '-q', '--bare', '--initial-branch=main', $remote);
command($dir, 'git', 'clone', '-q', $remote, $local);
command($local, 'git', 'config', 'user.name', 'Sync Test');
command($local, 'git', 'config', 'user.email', 'sync-test@example.invalid');
copy($script, "$local/update_github_upload.sh") or die "copy script: $!";
make_path("$local/MultiGWAS-Explorer");
write_file("$local/source.md", "Original source\n");
write_file("$local/MultiGWAS-Explorer/tracked_plot.png", "Original plot\n");
write_file("$local/.gitignore", "MultiGWAS-Explorer/cache/\n");
commit_all($local, 'Initial source and tracked result fixture');
command($local, 'git', 'push', '-q', 'origin', 'main');
command($dir, 'git', 'clone', '-q', $remote, $other);
command($other, 'git', 'config', 'user.name', 'Other Computer');
command($other, 'git', 'config', 'user.email', 'other-test@example.invalid');

remote_edit('remote.md', "Update from another computer\n");
my ($status, $output) = synchronize();
is($status, 0, 'default invocation downloads updates with no local edits') or diag $output;
is(read_file("$local/remote.md"), "Update from another computer\n", 'remote source arrived locally');
is(command($local, 'git', 'rev-parse', 'HEAD'), command($remote, 'git', 'rev-parse', 'main'),
    'clean local branch matches the remote');

# Preserve an existing user stash, staged edits, ignored data, and untracked data.
write_file("$local/source.md", "Older user stash\n");
command($local, 'git', 'stash', 'push', '-qm', 'Existing user stash');
my $user_stash = command($local, 'git', 'rev-parse', 'refs/stash');
write_file("$local/source.md", "Local staged edit\n");
command($local, 'git', 'add', 'source.md');
make_path("$local/MultiGWAS-Explorer/cache");
write_file("$local/MultiGWAS-Explorer/cache/gwas.gz", 'Keep ignored data');
write_file("$local/input.tsv", 'Keep untracked data');
remote_edit('remote.md', "Second remote update\n");
my $pull_head = command($remote, 'git', 'rev-parse', 'main');
($status, $output) = synchronize('--pull-only');
is($status, 0, 'pull-only integrates remote updates around staged local edits') or diag $output;
is(command($local, 'git', 'rev-parse', 'HEAD'), $pull_head, 'pull-only creates no local source commit');
is(command($remote, 'git', 'rev-parse', 'main'), $pull_head, 'pull-only does not push');
is(read_file("$local/source.md"), "Local staged edit\n", 'local edit restored after pull');
is(command($local, 'git', 'diff', '--cached', '--name-only'), 'source.md', 'staging state restored');
is(command($local, 'git', 'rev-parse', 'refs/stash'), $user_stash, 'preexisting user stash retained');
is(read_file("$local/MultiGWAS-Explorer/cache/gwas.gz"), 'Keep ignored data', 'ignored GWAS input retained');
is(read_file("$local/input.tsv"), 'Keep untracked data', 'untracked input retained');

# Local source edits plus remote commits must merge and upload; results stay local.
write_file("$local/new tool.pl", "print qq{new source};\n");
write_file("$local/MultiGWAS-Explorer/tracked_plot.png", "Changed local plot\n");
make_path("$local/MultiGWAS-Explorer/upload_manhattan_subset_test");
my $generated = 'MultiGWAS-Explorer/upload_manhattan_subset_test/intermediate.pl';
write_file("$local/$generated", "Generated runner\n");
command($local, 'git', 'add', $generated, 'MultiGWAS-Explorer/tracked_plot.png');
remote_edit('remote.md', "Third remote update\n");
($status, $output) = synchronize('Commit local source and synchronize');
is($status, 0, 'default mode merges and uploads concurrent source edits') or diag $output;
is(read_file("$local/remote.md"), "Third remote update\n", 'remote updates survive merge');
is(command($remote, 'git', 'show', 'main:source.md'), 'Local staged edit', 'local source uploaded');
is(command($remote, 'git', 'show', 'main:new tool.pl'), 'print qq{new source};', 'new source with spaces uploaded');
is(command($remote, 'git', 'show', 'main:MultiGWAS-Explorer/tracked_plot.png'), 'Original plot',
    'tracked generated result changes excluded');
is(read_file("$local/MultiGWAS-Explorer/tracked_plot.png"), "Changed local plot\n", 'excluded tracked result preserved locally');
is(read_file("$local/$generated"), "Generated runner\n", 'manually staged intermediate retained locally');
my ($missing_generated) = run_in($remote, 'git', 'cat-file', '-e', "main:$generated");
ok($missing_generated != 0, 'manually staged intermediate excluded from remote');
is(command($local, 'git', 'rev-parse', 'HEAD'), command($remote, 'git', 'rev-parse', 'main'),
    'merged history synchronized online');
is(command($local, 'git', 'rev-parse', 'refs/stash'), $user_stash, 'user stash retained after default synchronization');

# Ahead-only existing commits are uploaded even without new working tree edits.
write_file("$local/already_committed.md", "Already committed locally\n");
command($local, 'git', 'add', 'already_committed.md');
command($local, 'git', 'commit', '-qm', 'Existing local commit');
($status, $output) = synchronize();
is($status, 0, 'existing local-only commits are pushed without new source edits') or diag $output;
is(command($remote, 'git', 'show', 'main:already_committed.md'), 'Already committed locally', 'existing local commit uploaded');

# A preview cannot fetch, stage in the real index, commit, or modify local files.
write_file("$local/source.md", "Unstaged preview edit\n");
command($local, 'git', 'add', $generated);
my $index_before = read_file("$local/.git/index");
my $head_before = command($local, 'git', 'rev-parse', 'HEAD');
my $tracking_before = command($local, 'git', 'rev-parse', 'origin/main');
command($local, 'git', 'remote', 'set-url', 'origin', "$dir/unreachable.git");
($status, $output) = synchronize('--dry-run');
is($status, 0, 'dry-run works with an unreachable remote and does not fetch') or diag $output;
is(read_file("$local/.git/index"), $index_before, 'dry-run preserves the real Git index');
is(command($local, 'git', 'rev-parse', 'HEAD'), $head_before, 'dry-run preserves branch history');
is(command($local, 'git', 'rev-parse', 'origin/main'), $tracking_before, 'dry-run preserves remote tracking refs');
is(read_file("$local/source.md"), "Unstaged preview edit\n", 'dry-run preserves local file contents');
($status, $output) = synchronize('--dry-run', '--pull-only');
is($status, 0, 'pull-only dry-run also works offline') or diag $output;
is(read_file("$local/.git/index"), $index_before, 'pull-only dry-run preserves staging');
command($local, 'git', 'remote', 'set-url', 'origin', $remote);
command($local, 'git', 'restore', '--staged', $generated);
command($local, 'git', 'restore', 'source.md');

# Pulling a changed copy of the running script must use the original invocation.
command($other, 'git', 'pull', '-q', '--ff-only');
write_file("$other/update_github_upload.sh", read_file($script) . "\necho UNEXPECTED_SECOND_SCRIPT_EXECUTION >&2\n");
commit_all($other, 'Update the running synchronization script');
command($other, 'git', 'push', '-q', 'origin', 'main');
($status, $output) = synchronize('--pull-only');
is($status, 0, 'script can download an update to itself while running') or diag $output;
unlike($output, qr/UNEXPECTED_SECOND_SCRIPT_EXECUTION/, 'downloaded script tail is not executed during the active invocation');

# Merge conflicts stop before push and keep both histories for resolution.
write_file("$local/source.md", "Conflicting local source\n");
remote_edit('source.md', "Conflicting remote source\n");
my $conflict_head = command($remote, 'git', 'rev-parse', 'main');
($status, $output) = synchronize('Conflicting local source commit');
ok($status != 0, 'source merge conflict reports failure');
is(command($remote, 'git', 'rev-parse', 'main'), $conflict_head, 'merge conflict never pushes');
ok(length command($local, 'git', 'ls-files', '--unmerged'), 'conflicting files remain available for resolution');
is(command($local, 'git', 'show', 'refs/stash:MultiGWAS-Explorer/tracked_plot.png'), 'Changed local plot',
    'excluded result is safely retained in stash during a failed merge');
like($output, qr/git stash apply --index/, 'failed merge prints saved-edit recovery command');
my $failed_head = command($local, 'git', 'rev-parse', 'HEAD');
($status, $output) = synchronize();
ok($status != 0, 'rerun refuses to change a merge in progress');
is(command($local, 'git', 'rev-parse', 'HEAD'), $failed_head, 'in-progress merge guard preserves local commits');

# Clean up only the disposable fixture merge, then test conflict restoring edits.
command($local, 'git', 'merge', '--abort');
command($local, 'git', 'reset', '--hard', 'origin/main');
write_file("$local/source.md", "Uncommitted local source\n");
command($local, 'git', 'add', 'source.md');
remote_edit('source.md', "New remote source\n");
my $restore_conflict_head = command($remote, 'git', 'rev-parse', 'main');
($status, $output) = synchronize('--pull-only');
ok($status != 0, 'conflict restoring staged edits stops synchronization');
is(command($remote, 'git', 'rev-parse', 'main'), $restore_conflict_head, 'restore conflict never pushes');
is(command($local, 'git', 'show', 'refs/stash:source.md'), 'Uncommitted local source', 'conflicting local edits retained in stash');
like($output, qr/Your edits remain in stash/, 'stash conflict prints recovery information');

done_testing();
