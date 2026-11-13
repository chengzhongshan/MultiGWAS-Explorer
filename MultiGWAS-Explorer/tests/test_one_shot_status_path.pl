#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../DiffGWASDeps";
use SAS_ODA_Runner;
use File::Spec;
use File::Temp qw(tempdir);
use IPC::Open3;
use Symbol qw(gensym);
use Test::More;

# Execute the actual one-shot CLI with an intentionally failing SASPy stub.
# No SAS account, Java process or network connection is used.
my $root = tempdir('one shot status XXXX', TMPDIR => 1, CLEANUP => 1);
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
    my $pid = open3(my $in, my $out, $error, $^X,
        "$Bin/../run_sas_codes_or_script_in_ODA.pl", @args,
        '--output-prefix', File::Spec->catdir($root, $name), '--run-timeout-seconds', 20);
    close $in;
    my $captured = (do {local $/; <$out>} // '') . (do {local $/; <$error>} // '');
    waitpid($pid, 0);
    isnt($? >> 8, 0, "$name retains its failed exit status");
    like($captured, qr/POC_STUB_CONNECTION_FAILURE/, "$name reports the original connection failure");
    unlike($captured, qr/FileNotFoundError.*output\.run\.status/,
        "$name status path does not hide the connection failure");
}
done_testing();
