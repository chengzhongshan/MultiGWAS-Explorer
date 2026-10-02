#!/usr/bin/env perl
use strict;
use warnings;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);
use File::Path qw(remove_tree);
use Getopt::Long qw(GetOptions);
use JSON::PP qw(decode_json);

my ($workdir, $run_dir);
my (@outputs, @related_dirs);
GetOptions(
    'workdir=s'        => \$workdir,
    'run-dir=s'        => \$run_dir,
    'require-output=s@' => \@outputs,
    'related-dir=s@'  => \@related_dirs,
) or die "Usage: $0 --workdir DIR --run-dir DIR --require-output FILE [--related-dir DIR ...]\n";
die "--workdir, --run-dir, and at least one --require-output are required\n"
    unless defined($workdir) && defined($run_dir) && @outputs;

if (($ENV{KEEP_LOCAL_ODA_ARTIFACTS} // '') eq '1'
    || ($ENV{KEEP_RENDERED_DEBUG_FILES} // '') eq '1') {
    print "[cleanup] Keeping local SAS ODA helper folders for debugging.\n";
    exit 0;
}

my $root = abs_path($workdir) or die "Cannot resolve workdir $workdir\n";
my $run = direct_child($root, $run_dir);
my $run_name = basename($run);
die "Unsafe SAS ODA run-folder name: $run_name\n"
    unless $run_name =~ /^(?:run_manhattan_png|run_single_snp_with_gtf|run_local_hits_manhattan_png)_(\d{8}_\d{6})$/;
my $stamp = $1;
my $status_path = "$run/output.run.status.json";
if (!-s $status_path) {
    print "[cleanup] Keeping $run_name because its completed run status is unavailable.\n";
    exit 0;
}
open my $status_fh, '<', $status_path or die "Cannot read $status_path: $!\n";
local $/;
my $status = eval { decode_json(<$status_fh>) };
close $status_fh;
if (ref($status) ne 'HASH' || ($status->{state} // '') ne 'completed'
    || !$status->{success} || !$status->{complete}) {
    print "[cleanup] Keeping $run_name because SAS ODA did not report completed success.\n";
    exit 0;
}

for my $output (@outputs) {
    my $resolved = direct_child($root, $output);
    if (!-s $resolved) {
        print "[cleanup] Keeping $run_name because final output is empty: $resolved\n";
        exit 0;
    }
}

my @remove = ($run);
my %allowed = map { $_ => 1 } qw(
    upload_manhattan_png_macro
    upload_single_snp_with_gtf_support
    upload_local_hits_support
    .autogen_get_gtf_macro_local_mh
);
for my $candidate (@related_dirs) {
    next unless -e $candidate;
    my $path = direct_child($root, $candidate);
    my $name = basename($path);
    die "Unsafe related SAS ODA folder name: $name\n"
        unless $name =~ /^(.+)_(\d{8}_\d{6})$/
            && $allowed{$1} && $2 eq $stamp;
    push @remove, $path;
}

my $errors = [];
remove_tree(@remove, { safe => 1, error => \$errors });
if (@$errors) {
    my @messages;
    for my $entry (@$errors) {
        push @messages, map { "$_: $entry->{$_}" } keys %$entry;
    }
    die "Could not remove SAS ODA helper folders: " . join('; ', @messages) . "\n";
}
print "[cleanup] Removed successful SAS ODA helper folders: "
    . join(', ', map { basename($_) } @remove) . "\n";

sub direct_child {
    my ($root, $candidate) = @_;
    die "Refusing symlink path $candidate\n" if -l $candidate;
    my $resolved = abs_path($candidate)
        or die "Cannot resolve $candidate\n";
    die "Refusing path outside the SAS ODA workdir: $candidate\n"
        unless dirname($resolved) eq $root;
    return $resolved;
}
