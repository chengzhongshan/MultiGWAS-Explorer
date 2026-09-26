package FixedEffectMeta;
use strict;
use warnings;
use Exporter 'import';
use POSIX qw(erfc);

our @EXPORT_OK = qw(fixed_effect_meta);

# Both cohorts must contribute an allele-aligned effect and a positive SE.
# With nonzero rho, use the two-study generalized inverse-variance estimate.
sub fixed_effect_meta {
    my ($b1, $s1, $b2, $s2, $rho) = @_;
    $rho = 0 unless defined $rho;
    for my $value ($b1, $s1, $b2, $s2) {
        return unless defined $value && $value =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
    }
    return unless $s1 > 0 && $s2 > 0;
    die "Meta-analysis rho must be between -1 and 1\n"
        unless $rho =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/
        && $rho > -1 && $rho < 1;
    if ($rho != 0) {
        my $cov = $rho * $s1 * $s2;
        my $denom = $s1 * $s1 + $s2 * $s2 - 2 * $cov;
        return unless $denom > 0;
        my $beta = (($s2 * $s2 - $cov) * $b1
                  + ($s1 * $s1 - $cov) * $b2) / $denom;
        my $se = sqrt(($s1 * $s1 * $s2 * $s2 - $cov * $cov) / $denom);
        my $z = $beta / $se;
        return ($beta, $se, $z, erfc(abs($z) / sqrt(2)));
    }
    my $w1 = 1 / ($s1 * $s1);
    my $w2 = 1 / ($s2 * $s2);
    my $weight = $w1 + $w2;
    return unless $weight > 0;
    my $beta = ($b1 * $w1 + $b2 * $w2) / $weight;
    my $se = sqrt(1 / $weight);
    my $z = $beta / $se;
    my $p = erfc(abs($z) / sqrt(2));
    return ($beta, $se, $z, $p);
}

1;
