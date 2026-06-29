#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use JSON::PP qw(decode_json);

my $config = '';
GetOptions('config=s' => \$config) or die "Usage: perl emit_diff_gwas_runner_env.pl --config runner.json\n";
die "--config is required\n" unless length $config;
open my $fh, '<', $config or die "Cannot read config $config: $!\n";
local $/;
my $json = <$fh>;
close $fh;
my $cfg = decode_json($json);
die "Runner config root must be a JSON object\n" unless ref($cfg) eq 'HASH';

for my $key (sort keys %{$cfg}) {
    next unless $key =~ /^[A-Z][A-Z0-9_]*$/;
    my $value = $cfg->{$key};
    next unless defined $value;
    if (ref($value) eq 'ARRAY') {
        $value = join(' ', @{$value});
    }
    elsif (ref($value)) {
        next;
    }
    print "export $key=", shell_quote($value), "\n";
}

sub shell_quote {
    my ($text) = @_;
    $text = '' unless defined $text;
    $text =~ s/'/'"'"'/g;
    return "'$text'";
}
