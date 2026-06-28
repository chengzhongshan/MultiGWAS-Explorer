#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);
use lib $Bin;
use DiffGWASConfig qw(
  load_config_file
  normalize_pair_map
  ordered_prefixes
  default_base_cols
  default_value_fields
);

my $config_file = '';
my $dataset = 'scz_mh';
my $source_type = 'gzip';
my $remote_basename = '';
my $pair_col = 'PAIR_TAG';
my $base_cols = join(',', default_base_cols());
my $value_fields = join(',', default_value_fields());
my $pair_map = 'SCZ_W3_ALL_SEX=ALL,SCZ_W3_ASN_SEX=ASN,SCZ_W3_EUR_SEX=EUR';
my $prefix_order = 'ALL,ASN,EUR';

GetOptions(
    'config=s'          => \$config_file,
    'dataset=s'         => \$dataset,
    'source-type=s'     => \$source_type,
    'remote-basename=s' => \$remote_basename,
    'pair-col=s'        => \$pair_col,
    'base-cols=s'       => \$base_cols,
    'value-fields=s'    => \$value_fields,
    'pair-map=s'        => \$pair_map,
    'prefix-order=s'    => \$prefix_order,
) or die usage();

die "--remote-basename is required\n" unless length $remote_basename;
die "--source-type must be gzip or plain\n" unless $source_type eq 'gzip' || $source_type eq 'plain';

my $cfg = load_config_file($config_file);
$pair_col = cfg_or($cfg, 'pair_col', $pair_col);
$base_cols = cfg_or($cfg, 'base_cols', $base_cols);
$value_fields = cfg_or($cfg, 'value_fields', $value_fields);
$pair_map = $cfg->{pair_map} if exists $cfg->{pair_map};
$prefix_order = cfg_or($cfg, 'prefix_order', $prefix_order);
my $wide_columns = cfg_list($cfg, 'wide_columns');
my $char_lengths = cfg_hash($cfg, 'char_lengths');
my $alias_map = cfg_hash($cfg, 'alias_map');
my $post_alias_map = cfg_hash($cfg, 'post_alias_map');

my @base_cols = parse_list($base_cols);
my @value_fields = parse_list($value_fields);
my %pair_to_prefix = normalize_pair_map($pair_map);
my @prefix_order = parse_list($prefix_order);
my @pair_order = ordered_prefixes(\%pair_to_prefix, \@prefix_order);

my @out_cols;
if (@{$wide_columns}) {
    @out_cols = @{$wide_columns};
}
else {
    @out_cols = @base_cols;
    for my $prefix (@pair_order) {
        push @out_cols, map { "${prefix}_$_" } @value_fields;
    }
}

my %col_type = map { $_ => 'num' } @out_cols;
for my $col (keys %{$char_lengths}) {
    $col_type{$col} = 'char' if exists $col_type{$col};
}

my @length_specs;
my @input_specs;
for my $col (@out_cols) {
    if ($col_type{$col} eq 'char') {
        my $len = $char_lengths->{$col};
        push @length_specs, "$col \$$len";
        push @input_specs, "$col :\$$len.";
    }
    else {
        push @input_specs, $col;
    }
}

my %out_col_exists = map { $_ => 1 } @out_cols;
my @alias_lines;
for my $alias (sort keys %{$alias_map}) {
    my $source = $alias_map->{$alias};
    next unless defined $source && $out_col_exists{$source};
    push @alias_lines, "  $alias = $source;";
}

my @z_lines;
my %derived_z_exists;
for my $col (@out_cols) {
    next unless $col =~ /_BETA$/;
    (my $se_col = $col) =~ s/_BETA$/_SE/;
    next unless $out_col_exists{$se_col};
    (my $z_col = $col) =~ s/_BETA$/_Z/;
    $derived_z_exists{$z_col} = 1;
    push @z_lines, "  if $se_col>0 then $z_col = $col / $se_col;";
}

my @post_alias_lines;
for my $alias (sort keys %{$post_alias_map}) {
    my $source = $post_alias_map->{$alias};
    next unless defined $source;
    next unless $out_col_exists{$source} || $derived_z_exists{$source};
    push @post_alias_lines, "  $alias = $source;";
}

my $filename_stmt = $source_type eq 'gzip'
  ? qq{filename mhdata zip "~/$remote_basename" gzip;}
  : qq{filename mhdata "~/$remote_basename";};

print "$filename_stmt\n\n";
print "data $dataset;\n";
print "  infile mhdata dlm='09'x dsd firstobs=2 truncover lrecl=32767;\n";
print "  length ", join(' ', @length_specs), ";\n" if @length_specs;
print "  input\n";
for my $spec (@input_specs) {
    print "    $spec\n";
}
print "  ;\n";
print join("\n", @alias_lines), "\n" if @alias_lines;
print join("\n", @z_lines), "\n" if @z_lines;
print join("\n", @post_alias_lines), "\n" if @post_alias_lines;
print "run;\n";

sub cfg_or {
    my ($cfg, $key, $fallback) = @_;
    return $fallback unless exists $cfg->{$key};
    my $value = $cfg->{$key};
    return ref($value) eq 'ARRAY' ? join(',', @{$value}) : $value;
}

sub cfg_hash {
    my ($cfg, $key) = @_;
    my %out;
    if (exists $cfg->{$key} && ref($cfg->{$key}) eq 'HASH') {
        %out = %{ $cfg->{$key} };
    }
    return \%out;
}

sub cfg_list {
    my ($cfg, $key) = @_;
    return [] unless exists $cfg->{$key} && ref($cfg->{$key}) eq 'ARRAY';
    return $cfg->{$key};
}

sub parse_list {
    my ($text) = @_;
    return grep { length } map { s/^\s+|\s+$//gr } split /\s*,\s*/, ($text // '');
}

sub usage {
    return <<"USAGE";
Usage:
  perl generate_sas_wide_import_include.pl --config FILE.json --remote-basename file.tsv.gz [options]

Options:
  --dataset NAME          Output SAS dataset name. Default: scz_mh
  --source-type TYPE      gzip or plain. Default: gzip
  --remote-basename NAME  Uploaded remote data filename in SAS ODA
  --pair-col NAME         Optional override for pair column
  --base-cols LIST        Optional override for wide key columns
  --value-fields LIST     Optional override for value fields
  --pair-map LIST         Optional override for pair/prefix mapping
  --prefix-order LIST     Optional override for output prefix ordering

Config extras:
  wide_columns            Explicit ordered list of already-wide columns to read.
                          When present, the helper skips pair_map/value_fields
                          expansion and imports these columns exactly.
  char_lengths            JSON object of character column lengths, for example:
                          {"A1":8,"A2":8,"SNP":40}
  alias_map               JSON object of derived aliases, for example:
                          {"ALL_STD_P":"ALL_STD_DIFF_P"}
  post_alias_map          JSON object of derived aliases that should be applied
                          after auto-generated Z-score variables, for example:
                          {"ALL_STD_Z":"ALL_STD_DIFF_Z","ASN_Z":"ASN_EUR_GROUP1_Z"}
USAGE
}
