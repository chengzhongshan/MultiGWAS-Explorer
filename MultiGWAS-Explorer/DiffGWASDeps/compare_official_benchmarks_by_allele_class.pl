#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use List::Util qw(sum);

my ($expected, $gwama, $easystrata, $out);
GetOptions(
    'expected=s'   => \$expected,
    'gwama=s'      => \$gwama,
    'easystrata=s' => \$easystrata,
    'out=s'        => \$out,
) or die "Invalid options\n";
die "Usage: $0 --expected expected.tsv --gwama gwama.out --easystrata results.txt --out comparison.tsv\n"
    unless defined($expected) && defined($gwama) && defined($easystrata) && defined($out);

sub table_header {
    my ($fh, $file) = @_;
    my $header = <$fh> // die "Empty file: $file\n";
    chomp $header;
    $header =~ s/\r$//;
    my @names = split /\t/, $header, -1;
    my %idx;
    @idx{@names} = (0 .. $#names);
    return %idx;
}

sub read_p {
    my ($file, $id_name, $p_name) = @_;
    open my $fh, '<', $file or die "Cannot open $file: $!\n";
    my %idx = table_header($fh, $file);
    die "Missing required columns in $file\n" unless exists($idx{$id_name}) && exists($idx{$p_name});
    my %p;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        my @f = split /\t/, $line, -1;
        my ($id, $value) = @f[$idx{$id_name}, $idx{$p_name}];
        next unless defined($id) && length($id) && defined($value)
            && $value =~ /\A(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\z/;
        $p{$id} = 0 + $value;
    }
    close $fh;
    return \%p;
}

my $truth = read_p($expected, 'SNP', 'DIFF_P');
my $gwp = read_p($gwama, 'rs_number', 'q_p-value');

open my $efh, '<', $easystrata or die "Cannot open $easystrata: $!\n";
my %eidx = table_header($efh, $easystrata);
for my $required (qw(SNP A1 A2 EASYSTRATA_PDIFF)) {
    die "Missing $required in $easystrata\n" unless exists $eidx{$required};
}
my (%esp, %allele_class);
while (my $line = <$efh>) {
    chomp $line;
    $line =~ s/\r$//;
    my @f = split /\t/, $line, -1;
    my ($id, $a1, $a2, $p) = @f[@eidx{qw(SNP A1 A2 EASYSTRATA_PDIFF)}];
    next unless defined($id) && length($id) && defined($p)
        && $p =~ /\A(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\z/;
    $esp{$id} = 0 + $p;
    my $key = join '', sort(uc($a1 // ''), uc($a2 // ''));
    $allele_class{$id} = ($key eq 'AT' || $key eq 'CG')
        ? 'STRAND_AMBIGUOUS' : 'NON_AMBIGUOUS';
}
close $efh;

sub summarize {
    my ($software, $observed, $scope) = @_;
    my (@abs, @logabs);
    my %discordant = map { $_ => 0 } (0.05, 1e-5, 5e-8);
    my ($expected_rows, $matched) = (0, 0);
    for my $id (keys %$truth) {
        next if $scope ne 'ALL' && ($allele_class{$id} // '') ne $scope;
        ++$expected_rows;
        next unless exists $observed->{$id};
        ++$matched;
        my ($a, $b) = ($truth->{$id}, $observed->{$id});
        push @abs, abs($a - $b);
        push @logabs, abs((-log($a) / log(10)) - (-log($b) / log(10))) if $a > 0 && $b > 0;
        for my $threshold (keys %discordant) {
            ++$discordant{$threshold} if (($a < $threshold) xor ($b < $threshold));
        }
    }
    @abs = sort { $a <=> $b } @abs;
    @logabs = sort { $a <=> $b } @logabs;
    return join("\t",
        $software, $scope, $expected_rows, $matched, $expected_rows - $matched,
        @abs ? sum(@abs) / @abs : 'NA',
        @abs ? $abs[-1] : 'NA',
        @logabs ? $logabs[-1] : 'NA',
        @discordant{0.05, 1e-5, 5e-8}
    );
}

open my $ofh, '>', $out or die "Cannot write $out: $!\n";
print {$ofh} join("\t", qw(
    SOFTWARE ALLELE_CLASS EXPECTED_ROWS MATCHED_ROWS MISSING_ROWS
    MEAN_ABS_P_DIFF MAX_ABS_P_DIFF MAX_ABS_LOG10P_DIFF
    DISCORDANT_P_LT_0P05 DISCORDANT_P_LT_1E_5 DISCORDANT_P_LT_5E_8
)), "\n";
for my $scope (qw(ALL NON_AMBIGUOUS STRAND_AMBIGUOUS)) {
    print {$ofh} summarize('GWAMA_2.2.2_q_p-value', $gwp, $scope), "\n";
    print {$ofh} summarize('EasyStrata_8.6_CALCPDIFF', \%esp, $scope), "\n";
}
close $ofh;
print "Wrote $out\n";
