#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw($GunzipError);
use File::Path qw(make_path);

my ($input, $outdir);
my $candidate_p = 1e-5;
my $bin_width = 1e-4;
my $max_chisq = 10;
my $exclude_strand_ambiguous = 1;
my $max_eaf_abs_diff = 0.20;
GetOptions(
    'input=s'                    => \$input,
    'outdir=s'                   => \$outdir,
    'candidate-p=f'              => \$candidate_p,
    'exclude-strand-ambiguous!'  => \$exclude_strand_ambiguous,
    'max-eaf-abs-diff=f'         => \$max_eaf_abs_diff,
) or die "Invalid options\n";
die "Usage: $0 --input standardized_long.tsv.gz --outdir DIR [--candidate-p 1e-5] [--[no-]exclude-strand-ambiguous] [--max-eaf-abs-diff 0.2]\n"
    unless defined $input && defined $outdir;
make_path($outdir) unless -d $outdir;

my $fh;
if ($input =~ /\.gz\z/i) {
    $fh = IO::Uncompress::Gunzip->new($input, MultiStream => 1)
        or die "Cannot gunzip $input: $GunzipError\n";
} else {
    open $fh, '<', $input or die "Cannot open $input: $!\n";
}

my $header = <$fh> // die "Empty input: $input\n";
chomp $header;
$header =~ s/^#//;
$header =~ s/\r$//;
my @names = split /\t/, $header, -1;
my %idx;
@idx{@names} = (0 .. $#names);
for my $required (qw(CHR BP A1 A2 SNP PAIR_TAG GROUP1_GWAS_TAG GROUP2_GWAS_TAG GROUP1_Z GROUP2_Z DIFF_Z DIFF_P IS_CHRX)) {
    die "Missing required column $required\n" unless exists $idx{$required};
}
my $has_frequency_fields = !grep { !exists $idx{$_} }
    qw(GROUP1_FRQ_A GROUP1_FRQ_U GROUP2_FRQ_A GROUP2_FRQ_U);

open my $qfh, '>', "$outdir/unfiltered_genomic_inflation.tsv"
    or die "Cannot write QC output: $!\n";
print {$qfh} join("\t", qw(
    PAIR_TAG STATISTIC GWAS_TAG N_VALID_AUTOSOMAL N_INVALID_AUTOSOMAL
    MEDIAN_CHISQ_BIN_MIDPOINT LAMBDA_GC HISTOGRAM_BIN_WIDTH
    SOURCE_PAIR_ROWS ELIGIBLE_VARIANT_ROWS AUTOSOMAL_VARIANT_ROWS X_VARIANT_ROWS
    EXCLUDED_STRAND_AMBIGUOUS EXCLUDED_EAF_DIFF_GT_MAX MAX_EAF_ABS_DIFF INPUT_SCOPE
)), "\n";

open my $sfh, '>', "$outdir/unfiltered_difference_candidates_p_lt_1e_5.tsv"
    or die "Cannot write signal output: $!\n";
print {$sfh} join("\t", qw(PAIR_TAG CHR BP SNP A1 A2 DIFF_Z DIFF_P)), "\n";

my (%hist, %valid, %invalid, %counts, %group1_tag, %group2_tag);
my $line_number = 1;
my $nbins = int($max_chisq / $bin_width);
my $chi1_median = 0.4549364231195727;

sub median_from_histogram {
    my ($h, $n) = @_;
    return ('NA', 'NA') unless $n;
    my $target = int(($n + 1) / 2);
    my $cumulative = 0;
    for my $bin (0 .. $nbins) {
        $cumulative += $h->[$bin] // 0;
        if ($cumulative >= $target) {
            my $mid = ($bin + 0.5) * $bin_width;
            return ($mid, $mid / $chi1_median);
        }
    }
    die "Median exceeded histogram range; increase --max-chisq\n";
}

sub finish_pair {
    my ($pair) = @_;
    my @labels = ('GROUP1', 'GROUP2', 'DIFFERENCE');
    my @tags = ($group1_tag{$pair}, $group2_tag{$pair}, $pair);
    for my $j (0 .. 2) {
        my ($median, $lambda) = median_from_histogram($hist{$pair}[$j], $valid{$pair}[$j] // 0);
        print {$qfh} join("\t",
            $pair, $labels[$j], $tags[$j], $valid{$pair}[$j] // 0, $invalid{$pair}[$j] // 0,
            $median, $lambda, $bin_width, $counts{$pair}{source} // 0, $counts{$pair}{total} // 0,
            $counts{$pair}{autosomal} // 0, $counts{$pair}{x} // 0,
            $counts{$pair}{excluded_ambiguous} // 0, $counts{$pair}{excluded_frequency} // 0,
            $max_eaf_abs_diff, 'full_unfiltered_harmonization_eligible_autosomes_for_lambda'
        ), "\n";
    }
}

sub numeric {
    my ($x) = @_;
    return defined $x && $x =~ /\A[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\z/;
}

sub mean_numeric {
    my (@x) = grep { numeric($_) } @_;
    return undef unless @x;
    my $sum = 0;
    $sum += $_ for @x;
    return $sum / @x;
}

sub strand_ambiguous {
    my ($a1, $a2) = @_;
    my $key = join '', sort(uc($a1 // ''), uc($a2 // ''));
    return $key eq 'AT' || $key eq 'CG';
}

while (my $line = <$fh>) {
    ++$line_number;
    chomp $line;
    $line =~ s/\r$//;
    next unless length $line;
    my @f = split /\t/, $line, -1;
    my $pair = $f[$idx{PAIR_TAG}];
    $group1_tag{$pair} //= $f[$idx{GROUP1_GWAS_TAG}];
    $group2_tag{$pair} //= $f[$idx{GROUP2_GWAS_TAG}];
    ++$counts{$pair}{source};
    if ($exclude_strand_ambiguous && strand_ambiguous(@f[@idx{qw(A1 A2)}])) {
        ++$counts{$pair}{excluded_ambiguous};
        next;
    }
    if ($has_frequency_fields) {
        my $eaf1 = mean_numeric(@f[@idx{qw(GROUP1_FRQ_A GROUP1_FRQ_U)}]);
        my $eaf2 = mean_numeric(@f[@idx{qw(GROUP2_FRQ_A GROUP2_FRQ_U)}]);
        if (defined($eaf1) && defined($eaf2) && abs($eaf1 - $eaf2) > $max_eaf_abs_diff) {
            ++$counts{$pair}{excluded_frequency};
            next;
        }
    }
    ++$counts{$pair}{total};
    my $is_x = ($f[$idx{IS_CHRX}] // '') eq '1';
    if ($is_x) {
        ++$counts{$pair}{x};
    } else {
        ++$counts{$pair}{autosomal};
        my @z = @f[@idx{qw(GROUP1_Z GROUP2_Z DIFF_Z)}];
        for my $j (0 .. 2) {
            if (defined $z[$j] && $z[$j] =~ /\A[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\z/) {
                my $chisq = $z[$j] * $z[$j];
                my $bin = int($chisq / $bin_width);
                $bin = $nbins if $bin > $nbins;
                ++$hist{$pair}[$j][$bin];
                ++$valid{$pair}[$j];
            } else {
                ++$invalid{$pair}[$j];
            }
        }
    }
    my $p = $f[$idx{DIFF_P}];
    if (defined $p && $p =~ /\A(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?\z/ && $p < $candidate_p) {
        print {$sfh} join("\t", @f[@idx{qw(PAIR_TAG CHR BP SNP A1 A2 DIFF_Z DIFF_P)}]), "\n";
    }
    print STDERR "Processed $line_number rows\n" if $line_number % 2_000_000 == 0;
}
finish_pair($_) for sort keys %counts;
close $fh;
close $qfh;
close $sfh;
print "Wrote unfiltered QC outputs to $outdir\n";
