package StudyDesignFilters;
use strict;
use warnings;
use Exporter 'import';
our @EXPORT_OK = qw(cross_ancestry_filter_advisory);

sub cross_ancestry_filter_advisory {
    my ($spec) = @_;
    return '' unless ref($spec) eq 'HASH' && ref($spec->{pairs}) eq 'ARRAY';
    my @cross_pairs;
    for my $pair (@{ $spec->{pairs} }) {
        next unless ref($pair) eq 'HASH';
        my $first = _ancestries($pair->{group1_label});
        my $second = _ancestries($pair->{group2_label});
        next unless length($first) && length($second) && $first ne $second;
        push @cross_pairs, $pair->{label} || "$pair->{group1_label} vs $pair->{group2_label}";
    }
    return '' unless @cross_pairs;
    my $eaf = $spec->{max_eaf_abs_diff};
    $eaf = 'none' unless defined($eaf) && length($eaf);
    my $maf = $spec->{top_hit_maf_threshold};
    $maf = 0 unless defined($maf) && length($maf);
    my $strand = exists($spec->{exclude_strand_ambiguous})
        ? ($spec->{exclude_strand_ambiguous} ? 'on' : 'off') : 'on';
    return 'Cross-ancestry comparison (' . join('; ', @cross_pairs)
        . "): max_eaf_abs_diff=$eaf, top_hit_maf_threshold=$maf, "
        . "exclude_strand_ambiguous=$strand. Population-specific allele frequencies and low MAF may be real; "
        . 'positive EAF/MAF cutoffs can remove association signals. Check allele alignment, '
        . 'ancestry-matched frequency sources, and exclusion counts. A SNP absent from one input '
        . 'cannot produce a paired differential effect even with these filters off.';
}

sub _ancestries {
    my ($label) = @_;
    my %seen;
    for my $token (split /[^A-Z]+/, uc($label // '')) {
        $token = 'EAS' if $token eq 'ASN';
        $seen{$token} = 1 if $token =~ /^(?:AFR|AMR|EAS|EUR|LAT|SAS)$/;
    }
    return join('+', sort keys %seen);
}

1;
