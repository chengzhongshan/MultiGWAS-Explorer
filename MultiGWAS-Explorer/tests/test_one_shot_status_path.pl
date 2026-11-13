#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../DiffGWASDeps";
use SAS_ODA_Runner;
use File::Spec;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Copy qw(copy);
use Cwd qw(getcwd);
use IPC::Open3;
use Symbol qw(gensym);
use Test::More;

# Execute the actual one-shot CLI with an intentionally failing SASPy stub.
# No SAS account, Java process or network connection is used.
my $root = tempdir('one shot status XXXX', TMPDIR => 1, CLEANUP => 1);
# A configured checkout prepends its installed site-packages to PYTHONPATH.
# Copy the real CLI/modules into a fixture without that environment so the
# intended stub cannot be shadowed by real SASPy during an installer check.
my $fixture = File::Spec->catdir($root, 'runner');
my $deps = File::Spec->catdir($fixture, 'DiffGWASDeps');
make_path($deps);
my $cli = File::Spec->catfile($fixture, 'run_sas_codes_or_script_in_ODA.pl');
copy("$Bin/../run_sas_codes_or_script_in_ODA.pl", $cli) or die $!;
opendir my $modules, "$Bin/../DiffGWASDeps" or die $!;
for my $name (grep { /\.pm\z/ } readdir $modules) {
    copy("$Bin/../DiffGWASDeps/$name", File::Spec->catfile($deps, $name)) or die $!;
}
closedir $modules;
open my $stub, '>', File::Spec->catfile($root, 'saspy.py') or die $!;
print {$stub} "class SASsession:\n    def __init__(self, **kwargs):\n        raise RuntimeError('POC_STUB_CONNECTION_FAILURE')\n";
close $stub;
my $python = $ENV{PIPELINE_PYTHON_BIN} || 'python3';
local $ENV{PIPELINE_PYTHON_BIN} = $python;
local $ENV{PYTHONPATH} = SAS_ODA_Runner::_normalize_local_path_for_python($root, $python);
local $ENV{SASPY_CFGNAME} = 'linuxlocal';
local $ENV{SASPY_CFGFILE} = '';
local $ENV{SAS_ODA_AUTOLOAD_MACROS} = 0;
local $ENV{OPEN_RESULT} = 0;
for my $case (
    ['submit', '--code', '%put OFFLINE_STATUS_TEST;'],
    ['listing', '--dir4listing', '~'],
) {
    my ($name, @args) = @$case;
    my $error = gensym;
    my $original_dir = getcwd();
    chdir $fixture or die $!; # The installed source may be a read-only SIF.
    my $pid = open3(my $in, my $out, $error, $^X,
        $cli, @args,
        '--output-prefix', File::Spec->catdir($root, $name), '--run-timeout-seconds', 20);
    chdir $original_dir or die $!;
    close $in;
    my $captured = (do {local $/; <$out>} // '') . (do {local $/; <$error>} // '');
    waitpid($pid, 0);
    isnt($? >> 8, 0, "$name retains its failed exit status");
    like($captured, qr/POC_STUB_CONNECTION_FAILURE/, "$name reports the original connection failure");
    unlike($captured, qr/FileNotFoundError.*output\.run\.status/,
        "$name status path does not hide the connection failure");
}
done_testing();
