#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);
use IO::Compress::Gzip qw($GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);
use lib $Bin;
use DiffGWASConfig qw(
  load_config_file
  normalize_pair_map
  ordered_prefixes
  default_base_cols
  default_value_fields
);

my $config_file = '';
my $project_tag;
my $input;
my $target_snp = '';
my $window_bp;
my $output = '';
my $manifest = '';
my $output_dir;
my $pair_col;
my $base_cols;
my $value_fields;
my $pair_map;
my $prefix_order;
my $htsbin = "$Bin/../local/bin";
my $target_chr_hint;
my $target_bp_hint;

GetOptions(
    'config=s'      => \$config_file,
    'project-tag=s' => \$project_tag,
    'input=s'      => \$input,
    'target-snp=s' => \$target_snp,
    'window-bp=f'  => \$window_bp,
    'output=s'     => \$output,
    'manifest=s'   => \$manifest,
    'output-dir=s' => \$output_dir,
    'pair-col=s'   => \$pair_col,
    'base-cols=s'  => \$base_cols,
    'value-fields=s'=> \$value_fields,
    'pair-map=s'   => \$pair_map,
    'prefix-order=s'=> \$prefix_order,
    'htsbin=s'     => \$htsbin,
    'target-chr=s' => \$target_chr_hint,
    'target-bp=f'  => \$target_bp_hint,
) or die usage();

my $cfg = load_config_file($config_file);
$project_tag = pick_value($project_tag, cfg_or($cfg, 'project_tag', undef), 'PGC_SCZ');
$input = pick_value(
    $input,
    cfg_or($cfg, 'input', undef),
    '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.tsv.gz'
);
$target_snp = pick_value($target_snp, cfg_or($cfg, 'target_snp', undef), '');
$window_bp = pick_value($window_bp, cfg_or($cfg, 'window_bp', undef), 1e7);
$output = pick_value($output, cfg_or($cfg, 'single_output', undef), '');
$manifest = pick_value($manifest, cfg_or($cfg, 'single_manifest', undef), '');
$output_dir = pick_value($output_dir, cfg_or($cfg, 'output_dir', undef), '.');
$pair_col = pick_value($pair_col, cfg_or($cfg, 'pair_col', undef), 'PAIR_TAG');
$base_cols = pick_value($base_cols, cfg_or($cfg, 'base_cols', undef), join(',', default_base_cols()));
$value_fields = pick_value($value_fields, cfg_or($cfg, 'value_fields', undef), join(',', default_value_fields()));
$pair_map = pick_value($pair_map, (exists $cfg->{pair_map} ? $cfg->{pair_map} : undef), 'SCZ_W3_ALL_SEX=ALL,SCZ_W3_ASN_SEX=ASN,SCZ_W3_EUR_SEX=EUR');
$prefix_order = pick_value($prefix_order, cfg_or($cfg, 'prefix_order', undef), 'ALL,ASN,EUR');
$target_chr_hint = pick_value($target_chr_hint, cfg_or($cfg, 'target_chr', undef), undef);
$target_bp_hint = pick_value($target_bp_hint, cfg_or($cfg, 'target_bp', undef), undef);

die "--target-snp is required\n" unless length $target_snp;
die "Input file not found: $input\n" unless -s $input;
die "--window-bp must be non-negative\n" unless $window_bp >= 0;
die "Output directory not found: $output_dir\n" unless -d $output_dir;

my @base_cols = parse_list($base_cols);
my @value_fields = parse_list($value_fields);
die "At least one --base-cols entry is required\n" unless @base_cols;
die "At least one --value-fields entry is required\n" unless @value_fields;

my ($header, $idx_ref) = read_header($input);
my %idx = %{$idx_ref};

for my $required (@base_cols, $pair_col, @value_fields) {
    die "Required column $required not found in header\n" unless exists $idx{$required};
}

my %pair_to_prefix = normalize_pair_map($pair_map);
my @prefix_order = parse_list($prefix_order);
my @pair_order = ordered_prefixes(\%pair_to_prefix, \@prefix_order);
die "No pair tags found in --pair-map\n" unless @pair_order;

my ($raw_target_chr, $target_bp, $target_lookup_mode);
if (defined $target_chr_hint && $target_chr_hint ne '' && defined $target_bp_hint) {
    $raw_target_chr = $target_chr_hint;
    $target_bp = $target_bp_hint;
    $target_lookup_mode = 'hint';
}
else {
    ($raw_target_chr, $target_bp) = find_target($input, \%idx, $target_snp);
    $target_lookup_mode = 'scan';
}
die "Target SNP not found in input: $target_snp\n" unless defined $raw_target_chr && defined $target_bp;

my $target_chr = sas_chr($raw_target_chr);
my $safe_target = safe_name($target_snp);
my $safe_chr = safe_name($target_chr);
my $safe_bp = safe_name($target_bp);
my $safe_window = safe_name(sprintf('%.0f', $window_bp));
my $tag = safe_name($project_tag || 'diff_gwas');
my $stem = "${tag}_single_snp_${safe_target}_chr${safe_chr}_bp${safe_bp}_window_${safe_window}";

$output ||= "$output_dir/$stem.wide.tsv.gz";
$manifest ||= "$output_dir/$stem.manifest.tsv";

my @out_cols = @base_cols;
for my $prefix (@pair_order) {
    push @out_cols, map { "${prefix}_$_" } @value_fields;
}

my %stats = (
    rows_read             => 0,
    rows_in_window        => 0,
    groups_seen           => 0,
    groups_written        => 0,
    rows_missing_pair     => 0,
    duplicate_prefixes    => 0,
    unknown_pair_tags     => 0,
    target_chr            => $target_chr,
    target_bp             => $target_bp,
    target_snp            => $target_snp,
    window_bp             => sprintf('%.0f', $window_bp),
    output                => $output,
    manifest              => $manifest,
    target_lookup_mode    => $target_lookup_mode,
    region_query_mode     => 'stream',
);

my $in = open_window_reader($input, $raw_target_chr, $target_bp, $window_bp, \%stats);
my $out = open_writer($output);
print {$out} join("\t", @out_cols), "\n";

my $current_key = '';
my @bucket;
while (my $line = <$in>) {
    chomp $line;
    $line =~ s/\r$//;
    next if $line =~ /^\s*$/;
    $stats{rows_read}++;

    my @v = split /\t/, $line, -1;
    my $row_chr_raw = $v[ $idx{CHR} ] // '';
    my $row_chr = sas_chr($row_chr_raw);
    next unless defined $row_chr && $row_chr eq $target_chr;

    my $bp = numeric($v[ $idx{BP} ]);
    next unless defined $bp;
    next if $bp < ($target_bp - $window_bp) || $bp > ($target_bp + $window_bp);

    $stats{rows_in_window}++;
    my $key = join("\t", map { $v[ $idx{$_} ] // '' } @base_cols);
    if (@bucket && $key ne $current_key) {
        process_bucket(\@bucket, $out, \%idx, \%pair_to_prefix, \@pair_order, \@out_cols, \%stats);
        @bucket = ();
    }
    $current_key = $key;
    push @bucket, \@v;
}
process_bucket(\@bucket, $out, \%idx, \%pair_to_prefix, \@pair_order, \@out_cols, \%stats) if @bucket;

close $in  or die "Failed closing input $input: $!\n";
close $out or die "Failed closing output $output: $!\n";

open my $man, '>', $manifest or die "Cannot write $manifest: $!\n";
print {$man} join("\t", qw(METRIC VALUE)), "\n";
for my $metric (qw(
  input
  output
  manifest
  target_snp
  target_chr
  target_bp
  window_bp
  target_lookup_mode
  region_query_mode
  pair_col
  pair_prefixes
  pair_map
  rows_read
  rows_in_window
  groups_seen
  groups_written
  rows_missing_pair
  duplicate_prefixes
  unknown_pair_tags
)) {
    my $value =
        $metric eq 'input'         ? $input
      : $metric eq 'output'        ? $output
      : $metric eq 'manifest'      ? $manifest
      : $metric eq 'pair_col'      ? $pair_col
      : $metric eq 'pair_prefixes' ? join(',', @pair_order)
      : $metric eq 'pair_map'      ? join(',', map { "$_=$pair_to_prefix{$_}" } sort keys %pair_to_prefix)
      : $stats{$metric};
    print {$man} join("\t", $metric, $value), "\n";
}
print {$man} join("\t", 'columns', join(',', @out_cols)), "\n";
close $man or die "Failed closing manifest $manifest: $!\n";

print join("\n",
    "OUTPUT\t$output",
    "MANIFEST\t$manifest",
    "TARGET_SNP\t$target_snp",
    "TARGET_CHR\t$target_chr",
    "TARGET_BP\t$target_bp",
    "WINDOW_BP\t" . sprintf('%.0f', $window_bp),
    "TARGET_LOOKUP_MODE\t$stats{target_lookup_mode}",
    "REGION_QUERY_MODE\t$stats{region_query_mode}",
    "ROWS_IN_WINDOW\t$stats{rows_in_window}",
    "GROUPS_WRITTEN\t$stats{groups_written}",
), "\n";

sub find_target {
    my ($path, $idx, $target) = @_;
    my $fh = open_reader($path);
    <$fh>;
    my ($target_chr, $target_bp);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^\s*$/;
        my @v = split /\t/, $line, -1;
        my $snp = $v[ $idx->{SNP} ] // '';
        next unless $snp eq $target;
        $target_chr = $v[ $idx->{CHR} ];
        $target_bp = numeric($v[ $idx->{BP} ]);
        last if defined $target_chr && defined $target_bp;
    }
    close $fh;
    return ($target_chr, $target_bp);
}

sub process_bucket {
    my ($bucket, $out, $idx, $pair_to_prefix, $pair_order, $out_cols, $stats) = @_;
    return unless @{$bucket};

    $stats->{groups_seen}++;

    my %row = map { $_ => '' } @{$out_cols};
    for my $base (@base_cols) {
        my $value = $bucket->[0][ $idx->{$base} ] // '';
        $value = sas_chr($value) if $base eq 'CHR';
        $row{$base} = $value;
    }

    my %seen_prefix;
    for my $vals (@{$bucket}) {
        my $pair_tag = $vals->[ $idx->{$pair_col} ] // '';
        if ($pair_tag eq '') {
            $stats->{rows_missing_pair}++;
            next;
        }

        my $prefix = $pair_to_prefix->{$pair_tag};
        unless (defined $prefix && $prefix ne '') {
            $stats->{unknown_pair_tags}++;
            next;
        }
        $stats->{duplicate_prefixes}++ if $seen_prefix{$prefix}++;

        for my $field (@value_fields) {
            my $out_col = "${prefix}_$field";
            my $value = $vals->[ $idx->{$field} ] // '';
            my $numeric = numeric($value);
            $row{$out_col} = defined $numeric ? fmt($numeric) : '';
        }
    }

    print {$out} join("\t", map { defined $row{$_} ? $row{$_} : '' } @{$out_cols}), "\n";
    $stats->{groups_written}++;
}

sub pick_value {
    my ($cli_value, $config_value, $default_value) = @_;
    return $cli_value if defined $cli_value && $cli_value ne '';
    return $config_value if defined $config_value && $config_value ne '';
    return $default_value;
}

sub cfg_or {
    my ($cfg, $key, $fallback) = @_;
    return $fallback unless exists $cfg->{$key};
    my $value = $cfg->{$key};
    return ref($value) eq 'ARRAY' ? join(',', @{$value}) : $value;
}

sub parse_list {
    my ($text) = @_;
    return grep { length } map { s/^\s+|\s+$//gr } split /\s*,\s*/, ($text // '');
}

sub read_header {
    my ($path) = @_;
    my $fh = open_reader($path);
    my $h = <$fh>;
    close $fh;
    die "Input is empty: $path\n" unless defined $h;
    chomp $h;
    $h =~ s/\r$//;
    $h =~ s/^#//;
    my @cols = split /\t/, $h, -1;
    my %idx;
    for my $i (0 .. $#cols) {
        $idx{$cols[$i]} = $i;
    }
    return ($h, \%idx);
}

sub open_reader {
    my ($path) = @_;
    my $fh;
    if ($path =~ /\.gz$/i) {
        if (prefer_shell_gzip($path)) {
            my $cmd = 'gzip -dc ' . shell_quote($path);
            open $fh, '-|', $cmd or die "Cannot read gzip input $path with gzip -dc: $!\n";
        }
        else {
            $fh = IO::Uncompress::Gunzip->new($path)
              or die "Cannot open gzip input $path: $GunzipError\n";
        }
    }
    else {
        open $fh, '<', $path or die "Cannot read $path: $!\n";
    }
    return $fh;
}

sub open_window_reader {
    my ($path, $raw_chr, $bp, $window_bp, $stats) = @_;
    my $start = int($bp - $window_bp);
    $start = 1 if $start < 1;
    my $end = int($bp + $window_bp);

    my $tabix = resolve_hts_tool($htsbin, 'tabix');
    if (defined $tabix && has_tabix_index($path)) {
        my @regions = $stats->{target_lookup_mode} && $stats->{target_lookup_mode} eq 'hint'
          ? candidate_regions($raw_chr, $start, $end)
          : ("$raw_chr:$start-$end");
        my $cmd = join(' ', shell_quote($tabix), shell_quote($path), map { shell_quote($_) } @regions);
        open my $fh, '-|', $cmd or die "Cannot tabix-query $path: $!\n";
        $stats->{region_query_mode} = 'tabix';
        return $fh;
    }

    my $fh = open_reader($path);
    <$fh>;
    $stats->{region_query_mode} = 'stream';
    return $fh;
}

sub open_writer {
    my ($path) = @_;
    my $fh;
    if ($path =~ /\.gz$/i) {
        if (prefer_shell_gzip($path)) {
            my $cmd = 'gzip -c > ' . shell_quote($path);
            open $fh, '|-', $cmd or die "Cannot write gzip output $path with gzip -c: $!\n";
        }
        else {
            $fh = IO::Compress::Gzip->new($path)
              or die "Cannot open gzip output $path: $GzipError\n";
        }
    }
    else {
        open $fh, '>', $path or die "Cannot write $path: $!\n";
    }
    return $fh;
}

sub prefer_shell_gzip {
    my ($path) = @_;
    return 1 if defined $path && $path =~ m{^/mnt/};
    return 0;
}

sub resolve_hts_tool {
    my ($dir, $tool) = @_;
    for my $candidate (
        (defined $dir && length $dir ? ("$dir/$tool", "$dir/$tool.exe") : ()),
        $tool,
    ) {
        next unless defined $candidate && length $candidate;
        return $candidate if -x $candidate || command_exists($candidate);
    }
    return undef;
}

sub command_exists {
    my ($cmd) = @_;
    return 0 unless defined $cmd && length $cmd;
    return scalar(`command -v '$cmd' 2>/dev/null`) ? 1 : 0;
}

sub has_tabix_index {
    my ($path) = @_;
    return 1 if -e "$path.tbi" || -e "$path.csi";
    if ($path =~ /\.(?:gz|bgz|bgzip)$/i) {
        (my $stem = $path) =~ s/\.(?:gz|bgz|bgzip)$//i;
        return 1 if -e "$stem.tbi" || -e "$stem.csi";
    }
    return 0;
}

sub candidate_regions {
    my ($chr, $start, $end) = @_;
    my @candidates;
    my %seen;
    my @base = ($chr);
    if (defined $chr) {
        (my $no_chr = $chr) =~ s/^chr//i;
        push @base, $no_chr if defined $no_chr && $no_chr ne $chr;
        push @base, "chr$no_chr" if defined $no_chr && $no_chr ne '' && "chr$no_chr" ne $chr;
        if ($no_chr =~ /^(?:23|X)$/i) {
            push @base, 'X', '23', 'chrX', 'chr23';
        }
    }
    for my $base_chr (@base) {
        next unless defined $base_chr && $base_chr ne '';
        my $region = "$base_chr:$start-$end";
        next if $seen{$region}++;
        push @candidates, $region;
    }
    return @candidates;
}

sub shell_quote {
    my ($text) = @_;
    return "''" unless defined $text && length $text;
    $text =~ s/'/'\"'\"'/g;
    return "'$text'";
}

sub numeric {
    my ($x) = @_;
    return undef unless defined $x;
    return undef if $x eq '' || $x =~ /^(?:NA|NaN|null|\.)$/i;
    return undef unless $x =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
    return 0 + $x;
}

sub fmt {
    my ($x) = @_;
    return '' unless defined $x;
    return sprintf('%.10g', $x);
}

sub sas_chr {
    my ($chr) = @_;
    return '' unless defined $chr;
    $chr =~ s/^chr//i;
    return 23 if $chr =~ /^(?:X|23)$/i;
    return $chr;
}

sub safe_name {
    my ($text) = @_;
    $text = '' unless defined $text;
    $text =~ s/[^A-Za-z0-9._-]+/_/g;
    $text =~ s/^_+//;
    $text =~ s/_+$//;
    return length($text) ? $text : 'NA';
}

sub usage {
    return <<"USAGE";
Usage:
  perl extract_single_snp_wide_diff_gwas.pl --target-snp rs123 [options]

Options:
  --config FILE.json       Reusable comparison config
  --project-tag NAME       Prefix used for auto-generated output names
  --input FILE.tsv.gz       Standardized long differential GWAS table
  --target-snp RSID         Required target SNP identifier
  --window-bp FLOAT         +/- window around target SNP. Default: 1e7
  --output FILE.tsv.gz      Output wide local subset
  --manifest FILE.tsv       Output run summary
  --output-dir DIR          Directory for auto-generated output paths. Default: .
  --pair-col NAME           Pair/group column. Default: PAIR_TAG
  --base-cols LIST          Comma-delimited key columns. Default: CHR,BP,A1,A2,SNP
  --value-fields LIST       Comma-delimited measure columns to pivot
  --pair-map LIST           Comma-delimited PAIR_TAG=PREFIX mapping
  --prefix-order LIST       Optional comma-delimited preferred output prefix order
  --target-chr VALUE        Optional known target chromosome hint to avoid a full-file SNP scan
  --target-bp FLOAT         Optional known target BP hint to avoid a full-file SNP scan
  --htsbin DIR              Directory containing tabix. Default: ../local/bin

Behavior:
  - Finds the requested target SNP in the long standardized GWAS input.
  - Extracts all rows within the target chromosome and +/- window.
  - If a bgzip/tabix index is present, the window extraction step uses tabix instead of a second full-file scan.
  - Pivots the local region into the wide beta/SE/P layout expected by the SAS ODA GTF runner.
  - Writes a key-value manifest including TARGET_CHR and TARGET_BP.
  - JSON config values override the built-in defaults, and command-line flags can override the config.
  - This script intentionally ignores generic config keys like output/manifest from the
    genome-wide subset workflow; use single_output/single_manifest in JSON if you want
    config-driven explicit paths for the one-locus extraction.
USAGE
}
