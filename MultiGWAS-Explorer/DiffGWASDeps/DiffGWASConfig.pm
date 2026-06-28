package DiffGWASConfig;

use strict;
use warnings;
use Exporter qw(import);
use JSON::PP qw(decode_json);

our @EXPORT_OK = qw(
  load_config_file
  normalize_pair_map
  ordered_prefixes
  prefix_rank
  default_base_cols
  default_value_fields
  default_filter_fields
  default_char_lengths
);

sub default_base_cols {
    return qw(CHR BP A1 A2 SNP);
}

sub default_value_fields {
    return qw(
      FEMALE_BETA MALE_BETA DIFF_BETA
      FEMALE_SE   MALE_SE   DIFF_SE
      FEMALE_P    MALE_P    DIFF_P    STD_DIFF_Z    STD_DIFF_P
      FEMALE_FRQ_A FEMALE_FRQ_U MALE_FRQ_A MALE_FRQ_U
      FEMALE_INFO MALE_INFO
    );
}

sub default_filter_fields {
    return qw(FEMALE_P MALE_P DIFF_P STD_DIFF_P);
}

sub default_char_lengths {
    return (
      A1  => 8,
      A2  => 8,
      SNP => 40,
    );
}

sub load_config_file {
    my ($path) = @_;
    return {} unless defined $path && length $path;
    open my $fh, '<', $path or die "Cannot read config $path: $!\n";
    local $/;
    my $json = <$fh>;
    close $fh;
    my $cfg = decode_json($json);
    die "Config root must be a JSON object: $path\n" unless ref($cfg) eq 'HASH';
    return $cfg;
}

sub normalize_pair_map {
    my ($pair_map_input) = @_;
    my %pair_to_prefix;
    my %prefix_seen;

    if (!defined $pair_map_input) {
        return %pair_to_prefix;
    }

    if (ref($pair_map_input) eq 'HASH') {
        for my $pair_tag (sort keys %{$pair_map_input}) {
            my $prefix = $pair_map_input->{$pair_tag};
            die "Invalid pair_map prefix for $pair_tag\n" unless defined $prefix && length $prefix;
            die "Duplicate prefix in pair_map: $prefix\n" if $prefix_seen{$prefix}++;
            $pair_to_prefix{$pair_tag} = $prefix;
        }
        return %pair_to_prefix;
    }

    if (ref($pair_map_input)) {
        die "pair_map must be a JSON object or a comma-delimited string\n";
    }

    for my $entry (grep { length } split /\s*,\s*/, $pair_map_input) {
        my ($pair_tag, $prefix) = split /\s*=\s*/, $entry, 2;
        die "Invalid pair_map entry: $entry\n" unless defined $pair_tag && defined $prefix;
        die "Duplicate pair tag in pair_map: $pair_tag\n" if exists $pair_to_prefix{$pair_tag};
        die "Duplicate prefix in pair_map: $prefix\n" if $prefix_seen{$prefix}++;
        $pair_to_prefix{$pair_tag} = $prefix;
    }

    return %pair_to_prefix;
}

sub prefix_rank {
    my ($prefix, $prefix_order) = @_;
    if (defined $prefix_order && ref($prefix_order) eq 'ARRAY') {
        for my $i (0 .. $#{$prefix_order}) {
            return $i + 1 if defined $prefix_order->[$i] && $prefix eq $prefix_order->[$i];
        }
    }
    return 100000;
}

sub ordered_prefixes {
    my ($pair_to_prefix, $prefix_order) = @_;
    my %seen;
    return sort {
        prefix_rank($a, $prefix_order) <=> prefix_rank($b, $prefix_order)
          || $a cmp $b
    } grep { !$seen{$_}++ } values %{$pair_to_prefix};
}

1;
