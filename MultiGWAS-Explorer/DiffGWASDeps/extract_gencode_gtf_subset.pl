#!/usr/bin/env perl
use strict;
use warnings;

use FindBin qw($Bin);
use lib $Bin;
use File::Basename qw(basename);
use File::Path qw(make_path);
use File::Spec;
use Fcntl qw(:flock);
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
my $use_tabix = 1;
my $tabix_bin = $ENV{TABIX_BIN} // '';
my $bgzip_bin = $ENV{BGZIP_BIN} // '';

GetOptions(
    'reference-build=s'               => \$reference_build,
    'gtf-gz=s'                        => \$gtf_gz,
    'gtf-url=s'                       => \$gtf_url,
    'cache-dir=s'                     => \$cache_dir,
    'output=s'                        => \$output,
    'region=s@'                       => \@regions,
    'include-non-protein-coding!'     => \$include_non_protein_coding,
    'use-tabix!'                      => \$use_tabix,
    'tabix-bin=s'                     => \$tabix_bin,
    'bgzip-bin=s'                     => \$bgzip_bin,
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

my ($in, $using_tabix, $indexed_gtf, $tabix_query_count) =
  open_gtf_input(
      gtf_gz      => $gtf_gz,
      regions     => \@parsed_regions,
      use_tabix   => $use_tabix,
      tabix_bin   => $tabix_bin,
      bgzip_bin   => $bgzip_bin,
  );
if (!$in) {
    die "tabix extraction is required by default but an indexed GTF query could not be opened. "
      . "Install/provide bgzip and tabix, or explicitly request the slower compatibility path with --no-use-tabix.\n"
      if $use_tabix;
    $in = IO::Uncompress::Gunzip->new($gtf_gz)
      or die "Cannot open $gtf_gz via gunzip: $GunzipError\n";
    $using_tabix = 0;
}
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
    next if !$using_tabix && $done_chr{$chr_text};
    next unless exists $regions_by_chr{$chr_text};
    if (!$using_tabix && $start > $max_end_by_chr{$chr_text}) {
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
close $in or die "Failed while reading " . ($using_tabix ? "tabix query output" : $gtf_gz) . "\n";

print "OUTPUT\t$output\n";
print "ROWS\t$rows\n";
print "REFERENCE_BUILD\t$reference_build\n";
print "GTF_URL\t$gtf_url\n";
print "INCLUDE_NON_PROTEIN_CODING\t$include_non_protein_coding\n";
print "GTF_ACCESS_MODE\t" . ($using_tabix ? 'TABIX' : 'SEQUENTIAL_GZIP') . "\n";
print "GTF_INDEXED_PATH\t$indexed_gtf\n" if $using_tabix;
print "TABIX_QUERY_COUNT\t$tabix_query_count\n" if $using_tabix;

sub parse_region {
    my ($text) = @_;
    die "Invalid --region value: $text\n" unless defined $text && $text =~ /^([^:]+):(\d+):(\d+)$/;
    my ($chr, $start, $end) = ($1, $2, $3);
    ($start, $end) = ($end, $start) if $start > $end;
    my $chr_text = uc(normalize_chr_text($chr));
    $chr_text = 'X'  if $chr_text eq '23';
    $chr_text = 'Y'  if $chr_text eq '24';
    $chr_text = 'MT' if $chr_text eq '25';
    return {
        chr   => $chr_text,
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

sub open_gtf_input {
    my (%args) = @_;
    return unless $args{use_tabix};
    my $tabix = find_executable($args{tabix_bin}, 'tabix');
    my $bgzip = find_executable($args{bgzip_bin}, 'bgzip');
    unless ($tabix && $bgzip) {
        die "tabix and bgzip are required for GTF extraction. Put them in DiffGWASDeps, "
          . "the pipeline root, local/bin, or PATH; alternatively set --tabix-bin/--bgzip-bin.\n";
    }

    my $indexed = prepare_bgzf_gtf(
        source => $args{gtf_gz},
        tabix  => $tabix,
        bgzip  => $bgzip,
    );
    return unless $indexed && -s $indexed && -s "$indexed.tbi";

    open my $list_fh, '-|', $tabix, '-l', to_native_path($indexed)
      or do {
          die "[tabix] Could not list indexed GTF contigs from $indexed.\n";
      };
    my (@contigs, %contig_for_chr);
    while (my $line = <$list_fh>) {
        chomp $line;
        next unless length $line;
        push @contigs, $line;
        my $key = uc(normalize_chr_text($line));
        $key = 'X'  if $key eq '23';
        $key = 'Y'  if $key eq '24';
        $key = 'MT' if $key eq '25';
        $contig_for_chr{$key} //= $line;
    }
    close $list_fh or do {
        die "[tabix] Could not read indexed GTF contigs from $indexed.\n";
    };

    my @merged = merge_regions($args{regions});
    my @queries;
    for my $region (@merged) {
        my $contig = $contig_for_chr{$region->{chr}};
        unless (defined $contig) {
            die "[tabix] Indexed GTF has no contig matching chromosome $region->{chr}.\n";
        }
        push @queries, "$contig:$region->{start}-$region->{end}";
    }
    return unless @queries;

    open my $query_fh, '-|', $tabix, to_native_path($indexed), @queries
      or do {
          die "[tabix] Could not start indexed GTF query for $indexed.\n";
      };
    print STDERR "[tabix] Querying " . scalar(@queries) . " merged GTF interval(s) from $indexed\n";
    return ($query_fh, 1, $indexed, scalar(@queries));
}

sub merge_regions {
    my ($regions) = @_;
    my %by_chr;
    push @{ $by_chr{$_->{chr}} }, { %$_ } for @{$regions || []};
    my @merged;
    for my $chr (sort chromosome_sort keys %by_chr) {
        my @sorted = sort { $a->{start} <=> $b->{start} || $a->{end} <=> $b->{end} } @{ $by_chr{$chr} };
        for my $region (@sorted) {
            if (@merged && $merged[-1]{chr} eq $chr && $region->{start} <= $merged[-1]{end} + 1) {
                $merged[-1]{end} = $region->{end} if $region->{end} > $merged[-1]{end};
            }
            else {
                push @merged, { %$region };
            }
        }
    }
    return @merged;
}

sub chromosome_sort {
    my $an = chr_text_to_num($a);
    my $bn = chr_text_to_num($b);
    return defined($an) && defined($bn) ? $an <=> $bn : $a cmp $b;
}

sub prepare_bgzf_gtf {
    my (%args) = @_;
    my $source = $args{source};
    return $source if -s "$source.tbi";

    my $indexed = $source;
    $indexed =~ s/\.gz$/.bgz/i;
    $indexed .= '.bgz' if $indexed eq $source;
    print STDERR "[tabix] Preparing BGZF GTF cache and index: $indexed, compared it with $source\n";
    my $lock_path = "$indexed.lock";
    open my $lock, '>>', $lock_path
      or die "[tabix] Cannot open index lock $lock_path: $!\n";
    flock($lock, LOCK_EX)
      or die "[tabix] Cannot lock $lock_path: $!\n";

    my $source_mtime = (stat($source))[9] // 0;
    my $indexed_mtime = -s $indexed ? ((stat($indexed))[9] // 0) : 0;
    my $index_mtime = -s "$indexed.tbi" ? ((stat("$indexed.tbi"))[9] // 0) : 0;
    if (-s $indexed && -s "$indexed.tbi"
        && $indexed_mtime >= $source_mtime && $index_mtime >= $indexed_mtime) {
        close $lock;
        return $indexed;
    }

    my $tmp_bgz = "$indexed.tmp.$$";
    my $tmp_tbi = "$tmp_bgz.tbi";
    unlink $tmp_bgz if -e $tmp_bgz;
    unlink $tmp_tbi if -e $tmp_tbi;
    print STDERR "[tabix] Building one-time BGZF GTF cache and index: $indexed\n";

    # Some GTFs (notably GENCODE's "liftXX" builds) preserve the line order
    # of the genome build they were lifted FROM, so they are not sorted by
    # chromosome/position in the lifted-to coordinates. Tabix requires
    # sorted input, so decompress -> sort -> bgzip rather than decompress
    # -> bgzip directly.
    my $sort_bin = find_executable('', 'sort');
    die "Cannot find 'sort' on PATH; required to prepare GTF for tabix indexing.\n"
      unless $sort_bin;

    my $tmp_sorted = "$indexed.sorted.tmp.$$";
    unlink $tmp_sorted if -e $tmp_sorted;

    my $source_fh = IO::Uncompress::Gunzip->new($source)
      or die "Cannot open $source via gunzip while building tabix cache: $GunzipError\n";
    pipe(my $sort_reader, my $sort_writer) or die "Cannot create sort pipe: $!\n";
    my $sort_pid = fork();
    die "Cannot fork sort process: $!\n" unless defined $sort_pid;
    if ($sort_pid == 0) {
        close $sort_writer;
        open STDIN, '<&', $sort_reader or die "Cannot connect sort stdin: $!\n";
        exec { $sort_bin } $sort_bin, '-t', "\t", '-k1,1', '-k4,4n', '-o', $tmp_sorted;
        die "Cannot execute $sort_bin: $!\n";
    }
    close $sort_reader;
    my $buffer;
    while (1) {
        my $read = read($source_fh, $buffer, 1024 * 1024);
        die "Failed reading $source: $!\n" unless defined $read;
        last if $read == 0;
        print {$sort_writer} $buffer or die "Failed streaming GTF to sort: $!\n";
    }
    close $source_fh;
    close $sort_writer or die "Failed closing sort input pipe: $!\n";
    waitpid($sort_pid, 0);
    die "sort failed while preparing $tmp_sorted\n" unless $? == 0 && -s $tmp_sorted;

    open my $sorted_fh, '<:raw', $tmp_sorted
      or die "Cannot open sorted GTF $tmp_sorted: $!\n";
    pipe(my $reader, my $writer) or die "Cannot create bgzip pipe: $!\n";
    my $pid = fork();
    die "Cannot fork bgzip process: $!\n" unless defined $pid;
    if ($pid == 0) {
        close $writer;
        open STDIN, '<&', $reader or die "Cannot connect bgzip stdin: $!\n";
        open STDOUT, '>:raw', $tmp_bgz or die "Cannot create $tmp_bgz: $!\n";
        exec { $args{bgzip} } $args{bgzip}, '-c';
        die "Cannot execute $args{bgzip}: $!\n";
    }
    close $reader;
    while (1) {
        my $read = read($sorted_fh, $buffer, 1024 * 1024);
        die "Failed reading $tmp_sorted: $!\n" unless defined $read;
        last if $read == 0;
        print {$writer} $buffer or die "Failed streaming GTF to bgzip: $!\n";
    }
    close $sorted_fh;
    close $writer or die "Failed closing bgzip input pipe: $!\n";
    waitpid($pid, 0);
    unlink $tmp_sorted;
    die "bgzip failed while creating $tmp_bgz\n" unless $? == 0 && -s $tmp_bgz;

    system { $args{tabix} } $args{tabix}, '-f', '-p', 'gff', to_native_path($tmp_bgz);
    die "tabix failed while indexing $tmp_bgz\n" unless $? == 0 && -s $tmp_tbi;
    unlink $indexed if -e $indexed;
    unlink "$indexed.tbi" if -e "$indexed.tbi";
    rename $tmp_bgz, $indexed or die "Cannot install BGZF cache $indexed: $!\n";
    rename $tmp_tbi, "$indexed.tbi" or die "Cannot install tabix index $indexed.tbi: $!\n";
    close $lock;
    return $indexed;
}

sub find_executable {
    my ($explicit, $name) = @_;
    if (defined $explicit && length $explicit) {
        return $explicit if -f $explicit && -x $explicit;
        warn "[tabix] Configured $name executable is unavailable: $explicit\n";
        return;
    }
    my @suffixes = $^O =~ /MSWin32/i ? ('', '.exe', '.bat', '.cmd') : ('', '.exe');
    my @dirs = (
        $Bin,
        File::Spec->catdir($Bin, File::Spec->updir()),
        File::Spec->catdir($Bin, File::Spec->updir(), 'local', 'bin'),
        File::Spec->path(),
    );
    my %seen;
    for my $dir (grep { defined($_) && length($_) && !$seen{$_}++ } @dirs) {
        for my $suffix (@suffixes) {
            my $candidate = File::Spec->catfile($dir, $name . $suffix);
            return $candidate if -f $candidate && -x $candidate;
        }
    }
    return;
}

my $PATH_TRANSLATOR;  # cached path to wslpath/cygpath, computed once

sub to_native_path {
    # Native Windows tabix.exe (unlike a POSIX-aware bgzip that only
    # streams via stdin/stdout) receives an explicit file path and calls
    # into the Windows C runtime to open it, which does not understand
    # POSIX-style paths like /mnt/e/... or /cygdrive/e/.... Convert to a
    # Windows-style path (E:\...) before handing it to tabix.
    my ($path) = @_;
    return $path unless defined $path && length $path;
    return $path if $path =~ m{^[A-Za-z]:[\\/]};  # already Windows-style

    if (!defined $PATH_TRANSLATOR) {
        $PATH_TRANSLATOR = '';
        for my $bin (qw(wslpath cygpath)) {
            my $found = find_executable('', $bin);
            if ($found) {
                $PATH_TRANSLATOR = $found;
                last;
            }
        }
    }
    return $path unless length $PATH_TRANSLATOR;

    open my $translate_fh, '-|', $PATH_TRANSLATOR, '-w', $path
      or return $path;
    my $win_path = <$translate_fh> // '';
    close $translate_fh;
    chomp $win_path;
    return length $win_path ? $win_path : $path;
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
  --use-tabix / --no-use-tabix        Require/build a BGZF + tabix GTF index (default on).
                                      --no-use-tabix enables the slow sequential compatibility path.
  --tabix-bin FILE                    Override the tabix executable.
  --bgzip-bin FILE                    Override the bgzip executable.
USAGE
}
