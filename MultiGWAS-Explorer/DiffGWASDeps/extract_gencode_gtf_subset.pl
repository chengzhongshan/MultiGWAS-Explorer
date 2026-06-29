#!/usr/bin/env perl
use strict;
use warnings;

use FindBin qw($Bin);
use lib $Bin;
use File::Basename qw(basename);
use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);
use HTTP::Tiny;
use IO::Uncompress::Gunzip qw($GunzipError);
use GenomeBuildProfile qw(
  canonicalize_reference_build
  build_profile_for
);

my $gtf_gz = '';
my $gtf_url = '';
my $reference_build = 'hg38';
my $cache_dir = '';
my $output = '';
my @regions;
my $include_non_protein_coding = 0;

GetOptions(
    'reference-build=s'               => \$reference_build,
    'gtf-gz=s'                        => \$gtf_gz,
    'gtf-url=s'                       => \$gtf_url,
    'cache-dir=s'                     => \$cache_dir,
    'output=s'                        => \$output,
    'region=s@'                       => \@regions,
    'include-non-protein-coding!'     => \$include_non_protein_coding,
) or die usage();

die "--output is required\n" unless length $output;
die "At least one --region chr:start:end is required\n" unless @regions;
my $profile = build_profile_for(build => $reference_build)
  or die "Unsupported --reference-build value '$reference_build'. Use hg19, hg38, or t2t.\n";
$reference_build = $profile->{build};
$gtf_url ||= $profile->{gtf_url};

my @parsed_regions = map { parse_region($_) } @regions;
my (%regions_by_chr, %max_end_by_chr, %done_chr);
for my $region (@parsed_regions) {
    push @{ $regions_by_chr{ $region->{chr} } }, $region;
    my $end = $region->{end};
    $max_end_by_chr{ $region->{chr} } = $end
      if !exists $max_end_by_chr{ $region->{chr} } || $end > $max_end_by_chr{ $region->{chr} };
}
my $remaining_target_chrs = scalar keys %regions_by_chr;

if (!$gtf_gz) {
    $cache_dir ||= '.cache';
    make_path($cache_dir) unless -d $cache_dir;
    my $basename = basename($gtf_url);
    die "Could not infer basename from --gtf-url\n" unless length $basename;
    $gtf_gz = "$cache_dir/$basename";
    download_if_missing($gtf_url, $gtf_gz);
}

die "GTF gzip file not found: $gtf_gz\n" unless -f $gtf_gz;

my $in = IO::Uncompress::Gunzip->new($gtf_gz)
  or die "Cannot open $gtf_gz via gunzip: $GunzipError\n";
open my $out, '>', $output or die "Cannot write $output: $!\n";

print {$out} join("\t", qw(
  chr
  chr_text
  chr_raw
  seqname
  source
  feature
  type
  start
  end
  st
  en
  bp1
  bp2
  txStart
  txEnd
  score
  strand
  frame
  gene
  gene_name
  gene_id
  transcript_name
  transcript_id
  gene_type
  transcript_type
  exon_id
  exon_number
  ensembl
  genesymbol
  protein_coding
  original_protein_coding
)), "\n";

my $rows = 0;
while (my $line = <$in>) {
    next if $line =~ /^\s*#/;
    chomp $line;
    my @f = split /\t/, $line, 9;
    next unless @f >= 9;
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = @f;
    next unless defined $feature && defined $start && defined $end;
    next unless $feature =~ /^(?:gene|transcript|exon)$/i;

    my $chr_text = normalize_chr_text($seqname);
    next unless length $chr_text;
    next if $done_chr{$chr_text};
    next unless exists $regions_by_chr{$chr_text};
    if ($start > $max_end_by_chr{$chr_text}) {
        $done_chr{$chr_text} = 1;
        $remaining_target_chrs--;
        last if $remaining_target_chrs <= 0;
        next;
    }
    next unless overlaps_any_region($chr_text, $start, $end, $regions_by_chr{$chr_text});

    my %attr = parse_gtf_attributes($attribute // '');
    my $gene_id = $attr{gene_id} // '';
    my $gene_name = $attr{gene_name} // '';
    my $transcript_id = $attr{transcript_id} // '';
    my $transcript_name = $attr{transcript_name} // '';
    my $gene_type = $attr{gene_type} // '';
    my $transcript_type = $attr{transcript_type} // '';
    my $gene_biotype = $attr{gene_biotype} // '';
    my $transcript_biotype = $attr{transcript_biotype} // '';
    my $exon_id = $attr{exon_id} // '';
    my $exon_number = $attr{exon_number} // '';
    my $gene = first_nonempty($gene_name, $gene_id, $transcript_name, $transcript_id, lc($feature));
    my $genesymbol = first_nonempty($gene_name, $gene, $transcript_name, $gene_id, $transcript_id, lc($feature));
    $gene_type = first_nonempty($gene_type, $gene_biotype);
    $transcript_type = first_nonempty($transcript_type, $transcript_biotype);
    my $bio_text = lc join(' ', grep { defined && length } $gene_type, $transcript_type, $gene_biotype, $transcript_biotype, $attribute // '');
    my $original_protein_coding = infer_protein_coding_flag(
        bio_text      => $bio_text,
        transcript_id => $transcript_id,
    );
    next if !$include_non_protein_coding && !$original_protein_coding;
    my $protein_coding = $include_non_protein_coding ? 1 : $original_protein_coding;
    my $chr_num = chr_text_to_num($chr_text);
    next unless defined $chr_num;

    print {$out} join("\t",
        $chr_num,
        $chr_text,
        $seqname,
        $seqname,
        ($source // ''),
        lc($feature),
        lc($feature),
        $start,
        $end,
        $start,
        $end,
        $start,
        $end,
        $start,
        $end,
        ($score // ''),
        ($strand // ''),
        ($frame // ''),
        sanitize($gene),
        sanitize($gene_name),
        sanitize($gene_id),
        sanitize($transcript_name),
        sanitize($transcript_id),
        sanitize($gene_type),
        sanitize($transcript_type),
        sanitize($exon_id),
        sanitize($exon_number),
        sanitize($source || 'gencode'),
        sanitize($genesymbol),
        $protein_coding,
        $original_protein_coding,
    ), "\n";
    $rows++;
}

close $out or die "Cannot close $output: $!\n";
close $in;

print "OUTPUT\t$output\n";
print "ROWS\t$rows\n";
print "REFERENCE_BUILD\t$reference_build\n";
print "GTF_URL\t$gtf_url\n";
print "INCLUDE_NON_PROTEIN_CODING\t$include_non_protein_coding\n";

sub parse_region {
    my ($text) = @_;
    die "Invalid --region value: $text\n" unless defined $text && $text =~ /^([^:]+):(\d+):(\d+)$/;
    my ($chr, $start, $end) = ($1, $2, $3);
    ($start, $end) = ($end, $start) if $start > $end;
    return {
        chr   => uc(normalize_chr_text($chr)),
        start => $start,
        end   => $end,
    };
}

sub normalize_chr_text {
    my ($value) = @_;
    return '' unless defined $value;
    $value =~ s/^\s+|\s+$//g;
    $value =~ s/^chr//i;
    return uc $value;
}

sub chr_text_to_num {
    my ($chr_text) = @_;
    return 23 if $chr_text eq 'X';
    return 24 if $chr_text eq 'Y';
    return 25 if $chr_text eq 'M' || $chr_text eq 'MT';
    return $chr_text =~ /^\d+$/ ? int($chr_text) : undef;
}

sub overlaps_any_region {
    my ($chr_text, $start, $end, $regions) = @_;
    my $norm_chr = uc(normalize_chr_text($chr_text));
    for my $region (@{$regions}) {
        next unless $norm_chr eq $region->{chr};
        return 1 if $end >= $region->{start} && $start <= $region->{end};
    }
    return 0;
}

sub parse_gtf_attributes {
    my ($text) = @_;
    my %attr;
    while ($text =~ /([A-Za-z0-9_]+)\s+"([^"]*)"/g) {
        $attr{$1} = $2 unless exists $attr{$1};
    }
    return %attr;
}

sub first_nonempty {
    for my $value (@_) {
        next unless defined $value;
        return $value if length $value;
    }
    return '';
}

sub infer_protein_coding_flag {
    my (%args) = @_;
    my $bio_text = lc($args{bio_text} // '');
    return 1 if index($bio_text, 'protein_coding') >= 0;
    my $txid = $args{transcript_id} // '';
    return 1 if $txid =~ /^(?:NM_|XM_)/i;
    return 0;
}

sub sanitize {
    my ($value) = @_;
    $value = '' unless defined $value;
    $value =~ s/[\r\n\t]/ /g;
    return $value;
}

sub download_if_missing {
    my ($url, $path) = @_;
    return if -s $path;
    my $http = HTTP::Tiny->new(
        timeout => 600,
        verify_SSL => 1,
    );
    my $response = $http->mirror($url, $path);
    die "Failed to download $url: $response->{status} $response->{reason}\n"
      unless $response->{success} || $response->{status} == 304;
    die "Downloaded file is empty: $path\n" unless -s $path;
}

sub usage {
    return <<"USAGE";
Usage:
  perl extract_gencode_gtf_subset.pl --output subset.tsv --region 5:1:1000000 [options]

Options:
  --reference-build BUILD              Select a built-in GTF profile: hg19, hg38, or t2t.
  --gtf-gz FILE                       Use an existing local GTF.gz file.
  --gtf-url URL                       Download URL when --gtf-gz is absent.
  --cache-dir DIR                     Cache directory for downloaded GTF.gz.
  --output FILE                       Output TSV path.
  --region CHR:START:END              Region to keep. May be repeated.
  --include-non-protein-coding        Keep non-protein-coding features. Default off.
  --no-include-non-protein-coding     Restrict to protein-coding features.
USAGE
}
