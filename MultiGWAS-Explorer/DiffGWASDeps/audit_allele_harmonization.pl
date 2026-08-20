#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw($GunzipError);
use File::Path qw(make_path);

my ($input, $outdir);
GetOptions('input=s' => \$input, 'outdir=s' => \$outdir)
    or die "Invalid options\n";
die "Usage: $0 --input merged_long.sorted.coord.tsv.gz --outdir DIR\n"
    unless defined $input && defined $outdir;
make_path($outdir) unless -d $outdir;

my $fh;
if ($input =~ /\.gz\z/i) {
    $fh = IO::Uncompress::Gunzip->new($input, MultiStream => 1)
        or die "Cannot gunzip $input: $GunzipError\n";
} else {
    open $fh, '<', $input or die "Cannot open $input: $!\n";
}
my $header = <$fh> // die "Empty input\n";
chomp $header;
$header =~ s/^#//;
$header =~ s/\r$//;
my @names = split /\t/, $header, -1;
my %idx;
@idx{@names} = (0 .. $#names);
for my $required (qw(CHR BP A1 A2 SNP GWAS_TAG FRQ_A FRQ_U BETA SE)) {
    die "Missing required column $required\n" unless exists $idx{$required};
}

my @pairs = qw(ALL EUR ASN);
my %stats;
my %sample_count;
my $sample_limit = 100;
open my $efh, '>', "$outdir/harmonization_exceptions_sample.tsv"
    or die "Cannot write exception sample: $!\n";
print {$efh} join("\t", qw(PAIR CHR BP SNP CLASSIFICATION FEMALE_A1 FEMALE_A2 MALE_A1 MALE_A2 EAF_ABS_DIFF)), "\n";

sub numeric {
    my ($x) = @_;
    return defined $x && $x =~ /\A[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\z/;
}

sub eaf {
    my ($r) = @_;
    my @x = grep { numeric($_) } ($r->{FRQ_A}, $r->{FRQ_U});
    return undef unless @x;
    my $sum = 0;
    $sum += $_ for @x;
    return $sum / @x;
}

sub complement {
    my ($a) = @_;
    return undef unless defined $a && length($a) == 1;
    my %c = (A => 'T', T => 'A', C => 'G', G => 'C');
    return $c{uc $a};
}

sub ambiguous_pair {
    my ($a1, $a2) = @_;
    my $key = join('', sort(uc($a1 // ''), uc($a2 // '')));
    return $key eq 'AT' || $key eq 'CG';
}

sub relation {
    my ($f, $m) = @_;
    my ($fa1, $fa2, $ma1, $ma2) = map { uc($_ // '') }
        ($f->{A1}, $f->{A2}, $m->{A1}, $m->{A2});
    return 'DIRECT' if $fa1 eq $ma1 && $fa2 eq $ma2;
    return 'SWAPPED' if $fa1 eq $ma2 && $fa2 eq $ma1;
    my ($c1, $c2) = (complement($ma1), complement($ma2));
    if (defined $c1 && defined $c2) {
        return 'COMPLEMENT_DIRECT' if $fa1 eq $c1 && $fa2 eq $c2;
        return 'COMPLEMENT_SWAPPED' if $fa1 eq $c2 && $fa2 eq $c1;
    }
    return 'INCOMPATIBLE';
}

sub sample_exception {
    my ($pair, $chr, $bp, $snp, $class, $f, $m, $ediff) = @_;
    return if ($sample_count{"$pair\t$class"} // 0) >= $sample_limit;
    ++$sample_count{"$pair\t$class"};
    print {$efh} join("\t", $pair, $chr, $bp, $snp, $class,
        $f ? ($f->{A1}, $f->{A2}) : ('NA', 'NA'),
        $m ? ($m->{A1}, $m->{A2}) : ('NA', 'NA'),
        defined $ediff ? $ediff : 'NA'), "\n";
}

sub evaluate_snp {
    my ($chr, $bp, $snp, $by_tag) = @_;
    for my $pair (@pairs) {
        my $ftag = "${pair}_FEMALE";
        my $mtag = "${pair}_MALE";
        my $frows = $by_tag->{$ftag} // [];
        my $mrows = $by_tag->{$mtag} // [];
        next if !@$frows && !@$mrows;
        ++$stats{$pair}{VARIANT_KEYS_EVALUATED};
        ++$stats{$pair}{DUPLICATE_TAG_KEYS} if @$frows > 1 || @$mrows > 1;
        if (!@$frows) {
            ++$stats{$pair}{ONLY_MALE};
            sample_exception($pair, $chr, $bp, $snp, 'ONLY_MALE', undef, $mrows->[0], undef);
            next;
        }
        if (!@$mrows) {
            ++$stats{$pair}{ONLY_FEMALE};
            sample_exception($pair, $chr, $bp, $snp, 'ONLY_FEMALE', $frows->[0], undef, undef);
            next;
        }
        ++$stats{$pair}{BOTH_PRESENT};
        my ($best_f, $best_m, $best_relation);
        my %rank = (DIRECT => 1, SWAPPED => 2, COMPLEMENT_DIRECT => 3, COMPLEMENT_SWAPPED => 4, INCOMPATIBLE => 5);
        for my $f (@$frows) {
            for my $m (@$mrows) {
                my $rel = relation($f, $m);
                if (!defined $best_relation || $rank{$rel} < $rank{$best_relation}) {
                    ($best_f, $best_m, $best_relation) = ($f, $m, $rel);
                }
            }
        }
        ++$stats{$pair}{$best_relation};
        my $ambiguous = ambiguous_pair($best_f->{A1}, $best_f->{A2});
        ++$stats{$pair}{STRAND_AMBIGUOUS} if $ambiguous;
        my ($feaf, $meaf) = (eaf($best_f), eaf($best_m));
        if (defined($meaf) && ($best_relation eq 'SWAPPED' || $best_relation eq 'COMPLEMENT_SWAPPED')) {
            $meaf = 1 - $meaf;
        }
        my $ediff = defined($feaf) && defined($meaf) ? abs($feaf - $meaf) : undef;
        ++$stats{$pair}{EAF_DIFF_GT_0P1} if defined($ediff) && $ediff > 0.1;
        ++$stats{$pair}{EAF_DIFF_GT_0P2} if defined($ediff) && $ediff > 0.2;
        ++$stats{$pair}{CURRENT_EXACT_ALLELE_MATCH} if $best_relation eq 'DIRECT';
        my $compatible = $best_relation ne 'INCOMPATIBLE';
        my $freq_ok = !defined($ediff) || $ediff <= 0.2;
        ++$stats{$pair}{REVISED_HARMONIZATION_ELIGIBLE} if $compatible && !$ambiguous && $freq_ok;
        if ($best_relation ne 'DIRECT' || $ambiguous || !$freq_ok) {
            my @flags = ($best_relation);
            push @flags, 'STRAND_AMBIGUOUS' if $ambiguous;
            push @flags, 'EAF_DIFF_GT_0P2' unless $freq_ok;
            sample_exception($pair, $chr, $bp, $snp, join('+', @flags), $best_f, $best_m, $ediff);
        }
    }
}

sub process_coordinate {
    my ($chr, $bp, $rows) = @_;
    return unless @$rows;
    my %by_snp;
    for my $r (@$rows) {
        push @{$by_snp{$r->{SNP}}{$r->{GWAS_TAG}}}, $r;
    }
    evaluate_snp($chr, $bp, $_, $by_snp{$_}) for keys %by_snp;
}

my ($current_chr, $current_bp) = ('', '');
my @coordinate_rows;
my $line_number = 1;
while (my $line = <$fh>) {
    ++$line_number;
    chomp $line;
    $line =~ s/\r$//;
    next unless length $line;
    my @f = split /\t/, $line, -1;
    my ($chr, $bp) = @f[@idx{qw(CHR BP)}];
    if (@coordinate_rows && ($chr ne $current_chr || $bp ne $current_bp)) {
        process_coordinate($current_chr, $current_bp, \@coordinate_rows);
        @coordinate_rows = ();
    }
    ($current_chr, $current_bp) = ($chr, $bp);
    my %r;
    @r{qw(A1 A2 SNP GWAS_TAG FRQ_A FRQ_U BETA SE)} = @f[@idx{qw(A1 A2 SNP GWAS_TAG FRQ_A FRQ_U BETA SE)}];
    push @coordinate_rows, \%r;
    print STDERR "Processed $line_number rows\n" if $line_number % 5_000_000 == 0;
}
process_coordinate($current_chr, $current_bp, \@coordinate_rows);
close $fh;
close $efh;

my @metrics = qw(
    VARIANT_KEYS_EVALUATED BOTH_PRESENT ONLY_FEMALE ONLY_MALE
    DIRECT SWAPPED COMPLEMENT_DIRECT COMPLEMENT_SWAPPED INCOMPATIBLE
    STRAND_AMBIGUOUS EAF_DIFF_GT_0P1 EAF_DIFF_GT_0P2 DUPLICATE_TAG_KEYS
    CURRENT_EXACT_ALLELE_MATCH REVISED_HARMONIZATION_ELIGIBLE
);
open my $ofh, '>', "$outdir/harmonization_audit.tsv" or die "Cannot write audit: $!\n";
print {$ofh} join("\t", 'PAIR', @metrics), "\n";
for my $pair (@pairs) {
    print {$ofh} join("\t", $pair, map { $stats{$pair}{$_} // 0 } @metrics), "\n";
}
close $ofh;
print "Wrote harmonization audit to $outdir\n";
