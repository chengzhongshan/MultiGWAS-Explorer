#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw(sum max);

my ($expected, $gwama, $easystrata, $out);
GetOptions(
    'expected=s'   => \$expected,
    'gwama=s'      => \$gwama,
    'easystrata=s' => \$easystrata,
    'out=s'        => \$out,
) or die "Invalid options\n";
for my $arg ($expected, $gwama, $easystrata, $out) {
    die "Usage: $0 --expected expected.tsv --gwama gwama.out --easystrata results.txt --out comparison.tsv\n"
        unless defined $arg && length $arg;
}

sub read_table {
    my ($file, $id_col, $p_col) = @_;
    open my $fh, '<', $file or die "Cannot open $file: $!\n";
    my $header = <$fh> // die "Empty file: $file\n";
    chomp $header;
    $header =~ s/\r$//;
    my @names = split /\t/, $header, -1;
    my %idx;
    @idx{@names} = (0 .. $#names);
    die "Missing $id_col in $file\n" unless exists $idx{$id_col};
    die "Missing $p_col in $file\n"  unless exists $idx{$p_col};
    my %p;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        my ($id, $value) = @f[$idx{$id_col}, $idx{$p_col}];
        next unless defined $id && length $id && defined $value && $value =~ /\A(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\z/;
        die "Duplicate identifier $id in $file\n" if exists $p{$id};
        $p{$id} = 0 + $value;
    }
    close $fh;
    return \%p;
}

my $truth = read_table($expected,   'SNP',       'DIFF_P');
my $gwp   = read_table($gwama,      'rs_number', 'q_p-value');
my $esp   = read_table($easystrata, 'SNP',       'EASYSTRATA_PDIFF');

sub percentile {
    my ($sorted, $prob) = @_;
    return 'NA' unless @$sorted;
    my $i = int($prob * (@$sorted - 1) + 0.5);
    return $sorted->[$i];
}

sub summarize {
    my ($software, $observed) = @_;
    my (@abs, @logabs);
    my %discordant = map { $_ => 0 } (0.05, 1e-5, 5e-8);
    my $matched = 0;
    for my $id (keys %$truth) {
        next unless exists $observed->{$id};
        ++$matched;
        my ($a, $b) = ($truth->{$id}, $observed->{$id});
        push @abs, abs($a - $b);
        if ($a > 0 && $b > 0) {
            push @logabs, abs((-log($a) / log(10)) - (-log($b) / log(10)));
        }
        for my $t (keys %discordant) {
            ++$discordant{$t} if (($a < $t) xor ($b < $t));
        }
    }
    @abs = sort { $a <=> $b } @abs;
    @logabs = sort { $a <=> $b } @logabs;
    return {
        software => $software,
        expected_rows => scalar(keys %$truth),
        observed_rows => scalar(keys %$observed),
        matched_rows => $matched,
        missing_rows => scalar(keys %$truth) - $matched,
        mean_abs_p_diff => @abs ? sum(@abs) / @abs : 'NA',
        p95_abs_p_diff => percentile(\@abs, 0.95),
        max_abs_p_diff => @abs ? $abs[-1] : 'NA',
        max_abs_log10p_diff => @logabs ? $logabs[-1] : 'NA',
        discordant_p_lt_0p05 => $discordant{0.05},
        discordant_p_lt_1e_5 => $discordant{1e-5},
        discordant_p_lt_5e_8 => $discordant{5e-8},
    };
}

my @rows = (
    summarize('GWAMA_2.2.2_q_p-value', $gwp),
    summarize('EasyStrata_8.6_CALCPDIFF', $esp),
);
my @cols = qw(
    software expected_rows observed_rows matched_rows missing_rows
    mean_abs_p_diff p95_abs_p_diff max_abs_p_diff max_abs_log10p_diff
    discordant_p_lt_0p05 discordant_p_lt_1e_5 discordant_p_lt_5e_8
);
open my $oh, '>', $out or die "Cannot write $out: $!\n";
print {$oh} join("\t", @cols), "\n";
for my $row (@rows) {
    print {$oh} join("\t", map { $row->{$_} } @cols), "\n";
}
close $oh;
print "Wrote $out\n";
