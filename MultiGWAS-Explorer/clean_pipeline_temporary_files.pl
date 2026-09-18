#!/usr/bin/env perl
use strict;
use warnings;
use Cwd qw(abs_path);
use File::Spec;
use Getopt::Long qw(GetOptions);

my ($directory, $apply, $help, $hours) = ('.', 0, 0, 24);
GetOptions('workdir=s' => \$directory, 'apply' => \$apply,
    'min-age-hours=f' => \$hours, 'help' => \$help) or die "Use --help for usage.\n";
if ($help) {
    print <<'HELP';
Usage: perl clean_pipeline_temporary_files.pl [--workdir DIR] [--min-age-hours N] [--apply]

Preview recognized local temporary files (default: at least 24 hours old).
Add --apply to permanently delete the listed files. After ALL pipeline jobs
using this directory have stopped, --min-age-hours 0 includes recent files.
Run once per output/work directory. No SAS ODA remote files are touched.

Only known generated scripts, intermediate subsets, helper JSON files,
timestamped upload staging directories, and timestamped
run_local_hits_with_gtf directories are considered. Git-tracked files,
symlinks, final results outside those run directories, configs, inputs and
caches are preserved. Unknown files inside an upload directory preserve it.
HELP
    exit 0;
}
die "Invalid age or extra arguments\n" if $hours < 0 || @ARGV;
my $root = abs_path($directory) or die "Cannot resolve $directory\n";
die "Not a directory: $root\n" unless -d $root;
chdir $root or die "Cannot enter $root: $!\n";
my %tracked;
# Git is optional for standalone installations; inside a checkout, failure is fatal.
my $ancestor = $root;
my $in_git = 0;
while (1) {
    $in_git = 1 if -f File::Spec->catfile($ancestor, '.git')
        || -f File::Spec->catfile($ancestor, '.git', 'HEAD');
    my $parent = abs_path(File::Spec->catdir($ancestor, '..'));
    last if !defined($parent) || $parent eq $ancestor;
    $ancestor = $parent;
}
if ($in_git) {
    open my $git, '-|', 'git', 'ls-files', '-z', '--', '.' or die "Cannot list tracked files: $!\n";
    local $/ = "\0";
    while (my $name = <$git>) { chomp $name; $tracked{$name} = 1; }
    close $git or die "Cannot verify tracked files; refusing cleanup\n";
}
my $stamp = qr/\d{8}_\d{6}/;
sub temporary_file {
    my ($name) = @_;
    return $name =~ /\A(?:sas_inline_(?:code|runner|result)|sas_submit_(?:result|worker_result|worker_code))_[A-Za-z0-9]{4}\.(?:sas|py|json)\z/
        || $name =~ /\Asas_action_(?:result_[A-Za-z0-9]{4}\.json|runner_[A-Za-z0-9]{4}\.py)\z/
        || $name =~ /\A(?:auto_(?:gtf|wide)_import_(?:single_snp|local_hits_with_gtf)|run_sas_local_debug_local_top_hits_with_gtf|run_sas_oda_(?:local_top_hits_with_gtf|single_snp_with_gtf))\.$stamp(?:\.part\d+)?\.sas\z/
        || $name =~ /\A(?:local_gtf_subset_local_hits|target_snp_augmented_local_(?:gtf|mh))_$stamp\.tsv(?:\.gz)?\z/
        || $name =~ /\Asingle_snp_ld_augmented_rs\d+_$stamp\.tsv\.gz\z/;
}
my $cutoff = time - $hours * 3600;
sub eligible {
    my ($path) = @_;
    return 0 if -l $path || !-f $path || $tracked{$path};
    my @st = lstat $path;
    return @st && $st[9] <= $cutoff;
}
opendir my $dh, '.' or die "Cannot list $root: $!\n";
my @names = sort readdir $dh;
closedir $dh;
sub collect_run_tree {
    my ($top) = @_;
    my (@tree_files, @tree_dirs);
    my $valid = 1;
    my $walk;
    $walk = sub {
        my ($dir) = @_;
        my @dir_st = lstat $dir;
        if (!@dir_st || -l _ || !-d _ || $dir_st[9] > $cutoff) {
            $valid = 0;
            return;
        }
        opendir my $run_dh, $dir or die "Cannot list $dir: $!\n";
        my @children = grep { $_ ne '.' && $_ ne '..' } readdir $run_dh;
        closedir $run_dh;
        for my $child (@children) {
            my $path = "$dir/$child";
            if (-l $path) { $valid = 0; next; }
            if (-d $path) { $walk->($path); next; }
            if (-f $path && eligible($path)) { push @tree_files, $path; next; }
            $valid = 0;
        }
        push @tree_dirs, $dir;
    };
    $walk->($top);
    return $valid ? (\@tree_files, \@tree_dirs) : (undef, undef);
}

my (@files, @dirs);
for my $name (@names) {
    next if -l $name;
    if (temporary_file($name) && eligible($name)) { push @files, $name; next; }
    if (-d $name && $name =~ /\Arun_local_hits_with_gtf_$stamp(?:_part\d+)?\z/) {
        my ($run_files, $run_dirs) = collect_run_tree($name);
        if ($run_files) {
            push @files, @$run_files;
            push @dirs, @$run_dirs;
        }
        next;
    }
    next unless -d $name && $name =~ /\Aupload_(?:local_hits_with_gtf_manifest|top_hits_batches|manhattan_png_macro|manhattan_subset|forest_macros|forest_top_hits_csv|single_snp_with_gtf_support|local_hits_support|local_hits_subset|local_hits_requested_csv)_$stamp(?:_try\d+)?\z/;
    opendir my $stage, $name or die "Cannot list $name: $!\n";
    my @children = grep { $_ ne '.' && $_ ne '..' } readdir $stage;
    closedir $stage;
    my @paths = map { "$name/$_" } @children;
    # No recursive deletion: only a flat staging directory of known transfer files.
    next if grep { !/\.(?:sas|csv|tsv|tsv\.gz|json|txt|zip)\z/ || !eligible($_) } @paths;
    next if (stat($name))[9] > $cutoff;
    push @files, @paths;
    push @dirs, $name;
}
my $bytes = 0;
for my $path (sort @files) {
    die "File changed during cleanup: $path\n" unless eligible($path);
    $bytes += -s $path;
    print(($apply ? 'DELETE' : 'WOULD DELETE'), "\t$path\n");
    unlink $path or die "Cannot remove $path: $!\n" if $apply;
}
for my $path (@dirs) {
    print(($apply ? 'RMDIR' : 'WOULD RMDIR'), "\t$path\n");
    rmdir $path or die "Cannot remove empty $path: $!\n" if $apply;
}
printf "%s: %d files, %d directories, %.2f MiB.\n",
    ($apply ? 'Removed' : 'Preview'), scalar(@files), scalar(@dirs), $bytes / 1048576;
print "No files deleted. Use --apply after all pipeline jobs have stopped.\n" unless $apply;
