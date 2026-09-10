#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw(dirname);
use File::Spec;
use Cwd qw(abs_path);

# Generic entry point.  Keep the historical implementation filename so older
# commands remain valid while PLINK2 and HaploReg caches share one audit-tested
# clumping engine.
my $self = __FILE__;
if ($^O =~ /cygwin/i && $self =~ /^[A-Za-z]:[\\\/]/) {
    if (open my $cygpath, '-|', 'cygpath', '-u', $self) {
        my $converted = <$cygpath> // '';
        close $cygpath;
        $converted =~ s/[\r\n]+\z//;
        $self = $converted if length $converted;
    }
}
my $implementation = File::Spec->catfile(
    dirname(abs_path($self) || $self), 'clump_top_hits_with_haploreg_cache.pl'
);
exec $^X, $implementation, @ARGV;
die "Cannot run $implementation: $!\n";
