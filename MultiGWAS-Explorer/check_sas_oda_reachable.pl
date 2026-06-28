#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use File::Spec;
use POSIX qw(strftime);
use Cwd qw(abs_path getcwd);
use FindBin qw($Bin);

sub usage {
    return <<"USAGE";
Usage:
  perl check_sas_oda_reachable.pl [options]

Options:
  --output-prefix NAME        Output directory prefix for helper artifacts.
                              Default: sas_oda_reachability_YYYYMMDD_HHMMSS
  --timeout-seconds N         Helper submit timeout. Default: 60
  --grace-seconds N           Helper timeout grace period. Default: 10
  --session-id ID             Optional SAS ODA session id to probe with.
  --persistent                Reuse the provided session id if desired.
  --workdir DIR               Working directory for the helper call.
                              Default: current directory
  --help                      Show this help

This script is intentionally non-invasive: it submits only a tiny `%put`
statement through `run_sas_codes_or_script_in_ODA.pl` and reports whether SAS
ODA appears reachable from the current environment.
USAGE
}

my $output_prefix = 'sas_oda_reachability_' . strftime('%Y%m%d_%H%M%S', localtime);
my $timeout_seconds = 60;
my $grace_seconds = 10;
my $session_id = '';
my $persistent = 0;
my $workdir = getcwd();
my $help = 0;

GetOptions(
    'output-prefix=s'   => \$output_prefix,
    'timeout-seconds=i' => \$timeout_seconds,
    'grace-seconds=i'   => \$grace_seconds,
    'session-id=s'      => \$session_id,
    'persistent!'       => \$persistent,
    'workdir=s'         => \$workdir,
    'help!'             => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

my ($helper, $helper_uses_search) = resolve_oda_helper();
my $probe_code = q{%put ODA_REACHABLE_MARKER;};
my @cmd = ($^X);
push @cmd, '-S' if $helper_uses_search;
push @cmd, (
    $helper,
    '--code', $probe_code,
    '--output-prefix', $output_prefix,
);

if ($persistent || length $session_id) {
    push @cmd, '--persistent';
    push @cmd, '--session-id', ($session_id || 'mysession');
}

local $ENV{SAS_ODA_RUN_TIMEOUT_SECONDS} = $timeout_seconds;
local $ENV{SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS} = $grace_seconds;

my $orig_dir = getcwd();
chdir $workdir or die "Cannot chdir to $workdir: $!\n";

my $cmd_string = join ' ', map { shell_quote($_) } @cmd;
my $cmd_output = qx{$cmd_string 2>&1};
my $exit_code = $? >> 8;

my $info_path = File::Spec->catfile($output_prefix, 'output.html.info.txt');
my $info_text = (-f $info_path) ? slurp($info_path) : '';
my $combined = join "\n", grep { defined && length } ($cmd_output, $info_text);

my @hard_fail_patterns = (
    qr/could not log on to the server/i,
    qr/server configuration is invalid/i,
    qr/SAS submit timed out/i,
    qr/ended without a readable result payload/i,
    qr/connection refused/i,
    qr/failed to establish/i,
    qr/no usable log/i,
);

my $reachable = 0;
if ($combined =~ /ODA_REACHABLE_MARKER/) {
    $reachable = 1;
}
elsif ($combined =~ /SAS Connection established/i && $exit_code == 0) {
    $reachable = 1;
}

for my $pat (@hard_fail_patterns) {
    if ($combined =~ $pat) {
        $reachable = 0;
        last;
    }
}

print "STATUS: " . ($reachable ? 'REACHABLE' : 'UNREACHABLE') . "\n";
print "EXIT_CODE: $exit_code\n";
print "WORKDIR: " . abs_path($workdir) . "\n";
print "OUTPUT_PREFIX: $output_prefix\n";
print "INFO_FILE: " . File::Spec->catfile(abs_path($workdir), $info_path) . "\n";
if (length $session_id) {
    print "SESSION_ID: $session_id\n";
}
print "\n";
print $cmd_output if length $cmd_output;

chdir $orig_dir if defined $orig_dir;
exit($reachable ? 0 : 1);

sub slurp {
    my ($path) = @_;
    open my $fh, '<', $path or return '';
    local $/;
    my $text = <$fh>;
    close $fh;
    return defined $text ? $text : '';
}

sub shell_quote {
    my ($text) = @_;
    $text = '' unless defined $text;
    $text =~ s/'/'"'"'/g;
    return "'$text'";
}

sub resolve_oda_helper {
    my @candidates = (
        File::Spec->catfile($Bin, 'DiffGWASDeps', 'run_sas_codes_or_script_in_ODA.pl'),
        File::Spec->catfile($Bin, 'run_sas_codes_or_script_in_ODA.pl'),
    );
    for my $candidate (@candidates) {
        return ($candidate, 0) if defined $candidate && -f $candidate;
    }
    return ('run_sas_codes_or_script_in_ODA.pl', 1);
}
