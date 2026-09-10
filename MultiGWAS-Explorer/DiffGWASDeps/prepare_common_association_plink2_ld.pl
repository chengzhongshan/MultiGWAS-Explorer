#!/usr/bin/env perl
use strict;
use warnings;

use File::Basename qw(dirname);
use File::Spec;
use Getopt::Long qw(GetOptions);

my ($spec, $runner_config, $verify_output, $candidates, $leads, $audit,
    $cache, $status, $pfile, $bfile, $plink2, $populations,
    $population_rule, $haploreg_cache);
my $r2_threshold = 0.1;
my $window_kb = 1000;
my $haploreg_cache_min_r2 = 0.2;
my $signal_threshold = 5e-8;
my $top_p_thresholds = '';

GetOptions(
    'spec=s'                    => \$spec,
    'runner-config=s'           => \$runner_config,
    'verify-output=s'           => \$verify_output,
    'candidates=s'              => \$candidates,
    'leads=s'                   => \$leads,
    'audit=s'                   => \$audit,
    'cache=s'                   => \$cache,
    'status=s'                  => \$status,
    'pfile=s'                   => \$pfile,
    'bfile=s'                   => \$bfile,
    'plink2=s'                  => \$plink2,
    'populations=s'             => \$populations,
    'population-rule=s'         => \$population_rule,
    'r2-threshold=f'            => \$r2_threshold,
    'window-kb=f'               => \$window_kb,
    'signal-threshold=f'        => \$signal_threshold,
    'top-p-thresholds=s'        => \$top_p_thresholds,
    'haploreg-cache=s'          => \$haploreg_cache,
    'haploreg-cache-min-r2=f'   => \$haploreg_cache_min_r2,
) or die usage();

die usage() unless $spec && $verify_output && $candidates && $leads
    && $audit && $cache && $status;
$plink2 ||= $ENV{PLINK2} || 'plink2';
$populations ||= 'EUR EAS';
$population_rule ||= 'ANY';

my $bin = dirname(File::Spec->rel2abs(__FILE__));
my $verifier = File::Spec->catfile($bin, 'verify_common_association_loci.pl');
my $selector = File::Spec->catfile($bin, 'select_ld_pruned_top_hits.pl');

for my $path ($verify_output, $candidates, $leads, $audit, $cache, $status) {
    my $parent = dirname($path);
    mkdir $parent unless -d $parent;
}

my @verify = ($^X, $verifier, '--spec', $spec,
    '--output', $verify_output, '--candidates-out', $candidates);
push @verify, '--runner-config', $runner_config if $runner_config;
push @verify, '--top-p-thresholds', $top_p_thresholds if length $top_p_thresholds;
run(@verify);

my @select = ($^X, $selector,
    '--candidates', $candidates,
    '--output-leads', $leads,
    '--output-audit', $audit,
    '--output-cache', $cache,
    '--output-status', $status,
    '--plink2', $plink2,
    '--signal-column', 'common_assoc_p',
    '--signal-threshold', $signal_threshold,
    '--populations', $populations,
    '--population-rule', $population_rule,
    '--r2-threshold', $r2_threshold,
    '--window-kb', $window_kb,
);
push @select, '--pfile', $pfile if $pfile;
push @select, '--bfile', $bfile if !$pfile && $bfile;
if ($haploreg_cache) {
    push @select, '--haploreg-cache', $haploreg_cache,
        '--haploreg-cache-min-r2', $haploreg_cache_min_r2;
}
run(@select);

print "COMMON_ASSOCIATION_VERIFY\t$verify_output\n";
print "COMMON_ASSOCIATION_CANDIDATES\t$candidates\n";
print "TOP_HIT_LD_LEADS\t$leads\n";
print "TOP_HIT_LD_AUDIT\t$audit\n";
print "TOP_HIT_LD_CACHE\t$cache\n";
print "TOP_HIT_LD_REFERENCE_STATUS\t$status\n";

sub run {
    my (@command) = @_;
    system(@command);
    die "Command failed (exit $?): @command\n" if $? != 0;
}

sub usage {
    return <<'USAGE';
Usage: prepare_common_association_plink2_ld.pl --spec SPEC
  --verify-output FILE --candidates FILE --leads FILE --audit FILE
  --cache FILE --status FILE [--runner-config FILE]
  [--pfile PREFIX | --bfile PREFIX] [--plink2 EXE]
  [--populations "EUR EAS"] [--population-rule ANY|ALL]
  [--r2-threshold 0.1] [--window-kb 1000]
  [--haploreg-cache FILE --haploreg-cache-min-r2 0.2]

Runs the common-association verifier followed by PLINK2/1000 Genomes LD
selection. HaploReg4 is used only when the local genotype reference is absent.
USAGE
}
