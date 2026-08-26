#!/usr/bin/env perl
use strict;
use warnings;
use JSON::PP ();
use POSIX qw(strftime);
use File::Basename qw(dirname);
use File::Path qw(make_path);

my $status_file = '';
while (@ARGV) {
    my $arg = shift @ARGV;
    if ($arg eq '--status-file') {
        $status_file = shift(@ARGV) // '';
        next;
    }
    last if $arg eq '--';
    die "Unknown background-runner option: $arg\n";
}
die "--status-file is required\n" unless length $status_file;
die "A command is required after --\n" unless @ARGV;

my $started = strftime('%Y-%m-%dT%H:%M:%SZ', gmtime());
my $raw = system { $ARGV[0] } @ARGV;
my $exit_code = $raw == -1 ? 127 : ($raw >> 8);
my $signal = $raw == -1 ? 0 : ($raw & 127);
my $finished = strftime('%Y-%m-%dT%H:%M:%SZ', gmtime());
my $status = {
    state       => 'complete',
    exit_code   => $exit_code,
    signal      => $signal,
    started_at  => $started,
    finished_at => $finished,
};

my $dir = dirname($status_file);
make_path($dir) if length($dir) && $dir ne '.' && !-d $dir;
my $tmp = "$status_file.tmp.$$";
open(my $fh, '>:raw', $tmp) or die "Cannot write $tmp: $!\n";
print {$fh} JSON::PP->new->canonical(1)->encode($status), "\n";
close($fh) or die "Cannot close $tmp: $!\n";
rename($tmp, $status_file) or die "Cannot promote $tmp to $status_file: $!\n";
exit $exit_code;
