package ChromosomeBounds;
use strict;
use warnings;
use Exporter 'import';
our @EXPORT_OK = qw(chromosome_length locus_window);

# Primary chromosome lengths from UCSC hg38, hg19, and hs1 chromosome sizes:
# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
# https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
# https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes.txt
my %lengths = (
    hg38 => [qw(248956422 242193529 198295559 190214555 181538259
        170805979 159345973 145138636 138394717 133797422 135086622
        133275309 114364328 107043718 101991189 90338345 83257441
        80373285 58617616 64444167 46709983 50818468 156040895 57227415)],
    hg19 => [qw(249250621 243199373 198022430 191154276 180915260
        171115067 159138663 146364022 141213431 135534747 135006516
        133851895 115169878 107349540 102531392 90354753 81195210
        78077248 59128983 63025520 48129895 51304566 155270560 59373566)],
    hs1 => [qw(248387328 242696752 201105948 193574945 182045439
        172126628 160567428 146259331 150617247 134758134 135127769
        133324548 113566686 101161492 99753195 96330374 84276897
        80542538 61707364 66210255 45090682 51324926 154259566 62460029)],
);

sub chromosome_length {
    my ($build, $chr) = @_;
    $build = lc($build // '');
    $build = 'hg38' if $build eq 'grch38';
    $build = 'hg19' if $build eq 'grch37';
    $build = 'hs1' if $build eq 't2t' || $build eq 'chm13' || $build eq 't2t-chm13v2.0';
    return undef unless exists $lengths{$build};
    $chr = uc($chr // '');
    $chr =~ s/^CHR//;
    $chr = 23 if $chr eq 'X';
    $chr = 24 if $chr eq 'Y';
    return undef unless $chr =~ /^\d+$/ && $chr >= 1 && $chr <= 24;
    return $lengths{$build}[$chr - 1];
}

sub locus_window {
    my ($build, $chr, $bp, $half_window) = @_;
    my $start = $bp - $half_window;
    $start = 1 if $start < 1;
    my $end = $bp + $half_window;
    my $length = chromosome_length($build, $chr);
    $end = $length if defined($length) && $end > $length;
    die "Target coordinate $chr:$bp exceeds $build chromosome length $length\n"
        if defined($length) && $bp > $length;
    return ($start, $end);
}

1;
