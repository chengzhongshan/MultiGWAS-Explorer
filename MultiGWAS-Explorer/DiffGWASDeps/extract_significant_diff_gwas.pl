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
  default_filter_fields
);

my $config_file = '';
my $project_tag = 'PGC_SCZ';
my $input =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.tsv.gz';
my $output = '';
my $manifest = '';
my $threshold = 0.05;
my $pair_col = 'PAIR_TAG';
my $filter_fields = join(',', default_filter_fields());
my $value_fields = join(',', default_value_fields());
my $pair_map = 'SCZ_W3_ALL_SEX=ALL,SCZ_W3_ASN_SEX=ASN,SCZ_W3_EUR_SEX=EUR';
my $prefix_order = 'ALL,ASN,EUR';
my $base_cols = join(',', default_base_cols());

GetOptions(
    'config=s'        => \$config_file,
    'project-tag=s'   => \$project_tag,
    'input=s'         => \$input,
    'output=s'        => \$output,
    'manifest=s'      => \$manifest,
    'threshold=f'     => \$threshold,
    'pair-col=s'      => \$pair_col,
    'base-cols=s'     => \$base_cols,
    'value-fields=s'  => \$value_fields,
    'filter-fields=s' => \$filter_fields,
    'pair-map=s'      => \$pair_map,
    'prefix-order=s'  => \$prefix_order,
) or die usage();

my $cfg = load_config_file($config_file);
$project_tag = cfg_or($cfg, 'project_tag', $project_tag);
$input = cfg_or($cfg, 'input', $input);
$output = cfg_or($cfg, 'output', $output);
$manifest = cfg_or($cfg, 'manifest', $manifest);
$threshold = cfg_or($cfg, 'threshold', $threshold);
$pair_col = cfg_or($cfg, 'pair_col', $pair_col);
$base_cols = cfg_or($cfg, 'base_cols', $base_cols);
$value_fields = cfg_or($cfg, 'value_fields', $value_fields);
$filter_fields = cfg_or($cfg, 'filter_fields', $filter_fields);
$pair_map = $cfg->{pair_map} if exists $cfg->{pair_map};
$prefix_order = cfg_or($cfg, 'prefix_order', $prefix_order);

die "Input file not found: $input\n" unless -s $input;
die "threshold must be positive\n" unless $threshold > 0;

my @base_cols = parse_list($base_cols);
my @value_fields = parse_list($value_fields);
die "At least one --base-cols entry is required\n" unless @base_cols;
die "At least one --value-fields entry is required\n" unless @value_fields;

my @filter_fields = grep { length } map { uc $_ } split /\s*,\s*/, $filter_fields;
die "At least one --filter-fields entry is required\n" unless @filter_fields;
my %pair_to_prefix = normalize_pair_map($pair_map);
die "At least one --pair-map entry is required\n" unless %pair_to_prefix;
my @prefix_order = parse_list($prefix_order);

if (!length $output || !length $manifest) {
    my $stem = auto_output_stem($input, $project_tag, $threshold);
    $output ||= "$stem.tsv.gz";
    $manifest ||= "$stem.manifest.tsv";
}

my $output_tmp = temp_output_path($output);
my $manifest_tmp = temp_output_path($manifest);

my ($header, $idx_ref) = read_header($input);
my %idx = %{$idx_ref};

for my $required (@base_cols, $pair_col) {
    die "Required column $required not found in header\n" unless exists $idx{$required};
}
for my $field (@value_fields) {
    die "Required column $field not found in header\n" unless exists $idx{$field};
}
for my $field (@filter_fields) {
    die "Filter column $field not found in header\n" unless exists $idx{$field};
}

my @pair_order = ordered_prefixes(\%pair_to_prefix, \@prefix_order);
die "No pair tags found in --pair-map\n" unless @pair_order;

my @out_cols = @base_cols;
for my $prefix (@pair_order) {
    push @out_cols, map { "${prefix}_$_" } @value_fields;
}

my %stats = (
    rows_read          => 0,
    groups_seen        => 0,
    rows_written       => 0,
    groups_skipped     => 0,
    rows_missing_pair  => 0,
    bad_filter_values  => 0,
    duplicate_prefixes => 0,
    unknown_pair_tags  => 0,
);

my $in  = open_reader($input);
my $out = open_writer($output_tmp);
print {$out} join("\t", @out_cols), "\n";

<$in>;

my $current_key = '';
my @bucket;
while (my $line = <$in>) {
    chomp $line;
    $line =~ s/\r$//;
    next if $line =~ /^\s*$/;
    my @v = split /\t/, $line, -1;
    my $key = join("\t", map { $v[ $idx{$_} ] // '' } @base_cols);
    if (@bucket && $key ne $current_key) {
        process_bucket(\@bucket, $out, \%idx, \%pair_to_prefix, \@pair_order, \@out_cols, \@filter_fields, \%stats, $threshold);
        @bucket = ();
    }
    $current_key = $key;
    push @bucket, \@v;
    $stats{rows_read}++;
}
process_bucket(\@bucket, $out, \%idx, \%pair_to_prefix, \@pair_order, \@out_cols, \@filter_fields, \%stats, $threshold) if @bucket;

close $in  or die "Failed closing input $input: $!\n";
close $out or die "Failed closing output $output_tmp: $!\n";

open my $man, '>', $manifest_tmp or die "Cannot write $manifest_tmp: $!\n";
print {$man} join("\t", qw(METRIC VALUE)), "\n";
for my $metric (qw(
  input
  output
  threshold
  pair_col
  filter_fields
  pair_prefixes
  pair_map
  rows_read
  groups_seen
  rows_written
  groups_skipped
  rows_missing_pair
  bad_filter_values
  duplicate_prefixes
  unknown_pair_tags
)) {
    my $value =
        $metric eq 'input'         ? $input
      : $metric eq 'output'        ? $output
      : $metric eq 'threshold'     ? $threshold
      : $metric eq 'pair_col'      ? $pair_col
      : $metric eq 'filter_fields' ? join(',', @filter_fields)
      : $metric eq 'pair_prefixes' ? join(',', @pair_order)
      : $metric eq 'pair_map'      ? join(',', map { "$_=$pair_to_prefix{$_}" } sort keys %pair_to_prefix)
      : $stats{$metric};
    print {$man} join("\t", $metric, $value), "\n";
}
print {$man} join("\t", 'columns', join(',', @out_cols)), "\n";
close $man;

install_atomic_output($output_tmp, $output);
install_atomic_output($manifest_tmp, $manifest);

print "Input:         $input\n";
print "Output:        $output\n";
print "Manifest:      $manifest\n";
print "Threshold:     $threshold\n";
print "Filter fields: ", join(',', @filter_fields), "\n";
print "Pair prefixes: ", join(',', @pair_order), "\n";
print "Rows read:     $stats{rows_read}\n";
print "Rows written:  $stats{rows_written}\n";
print "Groups seen:   $stats{groups_seen}\n";
print "Groups skipped:$stats{groups_skipped}\n";

sub process_bucket {
    my ($bucket, $out, $idx, $pair_to_prefix, $pair_order, $out_cols, $filter_fields, $stats, $threshold) = @_;
    return unless @{$bucket};

    $stats->{groups_seen}++;

    my %row = map { $_ => '' } @{$out_cols};
    for my $base (@base_cols) {
        my $value = $bucket->[0][ $idx->{$base} ] // '';
        $value = sas_chr($value) if $base eq 'CHR';
        $row{$base} = $value;
    }

    my %seen_prefix;
    my $keep = 0;

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
            $row{$out_col} = defined $numeric ? fmt($numeric) : $value;
        }

        for my $field (@{$filter_fields}) {
            my $numeric = numeric($vals->[ $idx->{$field} ]);
            if (defined $numeric) {
                $keep = 1 if $numeric < $threshold;
            }
            else {
                $stats->{bad_filter_values}++;
            }
        }
    }

    if ($keep) {
        print {$out} join("\t", map { defined $row{$_} ? $row{$_} : '' } @{$out_cols}), "\n";
        $stats->{rows_written}++;
    }
    else {
        $stats->{groups_skipped}++;
    }
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

sub auto_output_stem {
    my ($input_path, $project, $thr) = @_;
    my $dir = '.';
    my $base = $input_path;
    if ($input_path =~ m{^(.*)/([^/]+)$}) {
        $dir = $1;
        $base = $2;
    }
    $base =~ s/\.tsv(?:\.gz)?$//i;
    my $tag = safe_name($project || $base);
    my $thr_tag = safe_name(sprintf('%.6g', $thr));
    return "$dir/${tag}.wide_subset_p_lt_${thr_tag}.final";
}

sub temp_output_path {
    my ($path) = @_;
    if ($path =~ /(.*)(\.(?:gz|bgz|bgzip))$/i) {
        return $1 . '.tmp.' . $$ . $2;
    }
    return $path . '.tmp.' . $$;
}

sub install_atomic_output {
    my ($tmp, $final) = @_;
    unlink $final if -e $final;
    rename $tmp, $final or die "Cannot move $tmp to $final: $!\n";
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
  perl extract_significant_diff_gwas.pl [options]

Options:
  --config FILE.json          Reusable comparison config
  --project-tag NAME          Prefix used for auto-generated output names
  --input FILE.tsv.gz          Standardized differential GWAS table
  --output FILE.tsv.gz         Wide-format subset output
  --manifest FILE.tsv          Run summary
  --threshold FLOAT            Keep SNP if any selected P field < threshold. Default: 0.05
  --pair-col NAME              Pair/group column. Default: PAIR_TAG
  --base-cols LIST             Comma-delimited key columns. Default: CHR,BP,A1,A2,SNP
  --value-fields LIST          Comma-delimited measure columns to pivot
  --filter-fields LIST         Comma-delimited P fields used for filtering.
                               Default: FEMALE_P,MALE_P,DIFF_P,STD_DIFF_P
  --pair-map LIST              Comma-delimited PAIR_TAG=PREFIX mapping.
                               Default: SCZ_W3_ALL_SEX=ALL,SCZ_W3_ASN_SEX=ASN,
                               SCZ_W3_EUR_SEX=EUR
  --prefix-order LIST          Optional comma-delimited preferred output prefix order

Output columns:
  <base-cols>
  <PAIR>_<value-field> for each configured pair prefix and measure field

Notes:
  - Converts chrX to 23 for SAS-friendly chromosome coding.
  - Reads/writes .gz directly, so it works without zcat/gzip shell helpers.
  - Uses explicit group-to-prefix mapping so the transform can stream in one pass on large GWAS files.
  - JSON config values override the built-in defaults, and command-line flags can override the config.
USAGE
}
