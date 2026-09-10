#!/usr/bin/env perl
use strict;
use warnings;

use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Spec;
use Getopt::Long qw(GetOptions);

my ($candidates, $output_leads, $output_audit, $output_cache, $output_status,
    $pfile, $bfile, $plink2, $haploreg_cache, $sqlite, $snp_column,
    $chr_column, $bp_column, $signal_column, $signal_columns, $populations);
my $signal_threshold = 1;
my $r2_threshold = 0.1;
my $window_kb = 1000;
my $threads = 4;
my $population_rule = 'ANY';
my $haploreg_cache_min_r2 = 0.2;
my $max_leads = 0;

$plink2 = $ENV{PLINK2} || 'plink2';
$snp_column = 'SNP';
$chr_column = 'CHR';
$bp_column = 'BP';
$signal_column = '';
$signal_columns = '';
$populations = 'EUR EAS';

GetOptions(
    'candidates=s'            => \$candidates,
    'output-leads=s'          => \$output_leads,
    'output-audit=s'          => \$output_audit,
    'output-cache=s'          => \$output_cache,
    'output-status=s'         => \$output_status,
    'pfile=s'                 => \$pfile,
    'bfile=s'                 => \$bfile,
    'plink2=s'                => \$plink2,
    'haploreg-cache=s'        => \$haploreg_cache,
    'sqlite=s'                => \$sqlite,
    'snp-column=s'            => \$snp_column,
    'chr-column=s'            => \$chr_column,
    'bp-column=s'             => \$bp_column,
    'signal-column=s'         => \$signal_column,
    'signal-columns=s'        => \$signal_columns,
    'signal-threshold=f'      => \$signal_threshold,
    'populations=s'           => \$populations,
    'population-rule=s'       => \$population_rule,
    'r2-threshold=f'          => \$r2_threshold,
    'window-kb=f'             => \$window_kb,
    'threads=i'               => \$threads,
    'haploreg-cache-min-r2=f' => \$haploreg_cache_min_r2,
    'max-leads=i'             => \$max_leads,
) or die usage();

die usage() unless defined($candidates) && defined($output_leads)
    && defined($output_audit) && (length($signal_column) xor length($signal_columns));
die "Specify --pfile or --bfile, not both\n" if defined($pfile) && defined($bfile);

my $bin_dir = script_dir();
my $builder = File::Spec->catfile($bin_dir, 'build_plink2_1kg_ld_cache.pl');
my $clumper = File::Spec->catfile($bin_dir, 'clump_top_hits_with_ld_cache.pl');
my $direct_clumper = File::Spec->catfile($bin_dir, 'clump_top_hits_with_plink2_1kg.pl');
$output_cache ||= "$output_audit.plink2_cache.tsv";
$output_status ||= "$output_audit.reference_status.tsv";
$sqlite ||= "$output_audit.sqlite";

my $reference_available = defined($pfile) ? pfile_available($pfile)
                        : defined($bfile) ? bfile_available($bfile) : 0;
$reference_available = 0 unless executable_available($plink2);
my ($source, $cache_min, $cache_path, $fallback_bp);
if ($reference_available) {
    my @direct = ($^X, $direct_clumper, '--candidates', $candidates,
        '--plink2', $plink2, '--populations', $populations,
        '--population-rule', $population_rule, '--snp-column', $snp_column,
        '--chr-column', $chr_column, '--bp-column', $bp_column,
        '--signal-threshold', $signal_threshold, '--r2-threshold', $r2_threshold,
        '--window-kb', $window_kb, '--threads', $threads, '--max-leads', $max_leads,
        '--output-leads', $output_leads, '--output-audit', $output_audit,
        '--output-cache', $output_cache, '--output-status', $output_status);
    push @direct, '--signal-column', $signal_column if length $signal_column;
    push @direct, '--signal-columns', $signal_columns if length $signal_columns;
    push @direct, defined($pfile) ? ('--pfile', $pfile) : ('--bfile', $bfile);
    run(@direct);
    print "TOP_HIT_LD_SOURCE\tPLINK2_1KG_DIRECT\n";
    print "TOP_HIT_LD_REFERENCE_AVAILABLE\t1\nTOP_HIT_LD_STATUS\tCOMPLETE\n";
    exit 0;
} elsif (defined($haploreg_cache) && -s $haploreg_cache) {
    die "Requested r2 threshold $r2_threshold is below HaploReg cache coverage $haploreg_cache_min_r2\n"
        if $r2_threshold < $haploreg_cache_min_r2;
    ($source, $cache_min, $cache_path, $fallback_bp) =
        ('HAPLOREG4_FALLBACK', $haploreg_cache_min_r2, $haploreg_cache, 0);
    open my $status, '>:raw', $output_status
        or die "Cannot write $output_status: $!\n";
    print {$status} "status\tsource\treason\nFALLBACK\tHAPLOREG4\tPLINK2_1KG_REFERENCE_UNAVAILABLE\n";
    close $status;
} else {
    print STDERR "ERROR[20]: PLINK2 1KG reference is unavailable and no HaploReg fallback cache was supplied\n";
    exit 20;
}

my $clump_populations = $populations;
$clump_populations =~ s/\bEAS\b/ASN/ig if $source eq 'HAPLOREG4_FALLBACK';
my @clump = ($^X, $clumper, '--candidates', $candidates,
    '--cache', $cache_path, '--sqlite', $sqlite,
    '--output-leads', $output_leads, '--output-audit', $output_audit,
    '--snp-column', $snp_column, '--chr-column', $chr_column,
    '--bp-column', $bp_column, '--signal-threshold', $signal_threshold,
    '--populations', $clump_populations, '--population-rule', $population_rule,
    '--r2-threshold', $r2_threshold, '--cache-min-r2', $cache_min,
    '--cache-source', $source, '--fallback-distance-bp', $fallback_bp,
    '--max-leads', $max_leads);
push @clump, '--signal-column', $signal_column if length $signal_column;
push @clump, '--signal-columns', $signal_columns if length $signal_columns;
run(@clump);

print "TOP_HIT_LD_SOURCE\t$source\n";
print "TOP_HIT_LD_REFERENCE_AVAILABLE\t", ($reference_available ? 1 : 0), "\n";
print "TOP_HIT_LD_STATUS\tCOMPLETE\n";

sub pfile_available {
    my ($prefix) = @_;
    return -s "$prefix.pgen" && (-s "$prefix.pvar" || -s "$prefix.pvar.zst")
        && -s "$prefix.psam";
}

sub bfile_available {
    my ($prefix) = @_;
    return -s "$prefix.bed" && -s "$prefix.bim" && -s "$prefix.fam";
}

sub executable_available {
    my ($exe) = @_;
    return 1 if -f $exe;
    for my $dir (File::Spec->path()) {
        return 1 if -f File::Spec->catfile($dir, $exe);
        return 1 if $^O =~ /MSWin32|cygwin/i && -f File::Spec->catfile($dir, "$exe.exe");
    }
    return 0;
}

sub run {
    my (@cmd) = @_;
    system(@cmd);
    die "Command failed (exit $?): @cmd\n" if $? != 0;
}

sub script_dir {
    my $self = __FILE__;
    if ($^O =~ /cygwin/i && $self =~ /^[A-Za-z]:[\\\/]/) {
        if (open my $fh, '-|', 'cygpath', '-u', $self) {
            my $converted = <$fh> // '';
            close $fh;
            $converted =~ s/[\r\n]+\z//;
            $self = $converted if length $converted;
        }
    }
    return dirname(abs_path($self) || $self);
}

sub usage {
    return <<'USAGE';
Usage: select_ld_pruned_top_hits.pl --candidates FILE --output-leads FILE
  --output-audit FILE (--signal-column P | --signal-columns "P1 P2")
  [--pfile 1KG_PREFIX | --bfile 1KG_PREFIX] [--plink2 EXE]
  [--populations "EUR EAS"] [--r2-threshold 0.1] [--window-kb 1000]
  [--haploreg-cache FILE --haploreg-cache-min-r2 0.2]

PLINK2/1000 Genomes is always preferred when its complete reference triplet is
available. HaploReg is used only when that genotype reference is unavailable.
With the direct reference, non-estimable variants are retained and explicitly
reported; they are never converted to r2=0 or distance-pruned.
USAGE
}
