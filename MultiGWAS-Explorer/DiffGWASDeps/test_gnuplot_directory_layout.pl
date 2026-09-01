#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;

my $project_root = File::Spec->catdir($Bin, File::Spec->updir());
my $canonical_dir = File::Spec->catdir($Bin, 'gnuplot');
my $legacy_dir = File::Spec->catdir($Bin, 'gunplot');

die "Legacy implementation directory still exists: $legacy_dir\n"
    if -e $legacy_dir;
die "Canonical gnuplot implementation directory is missing: $canonical_dir\n"
    unless -d $canonical_dir;

for my $file (qw(
    pdl_gunplot_manhattan.pl
    pdl_gunplot_forest.pl
    pdl_gunplot_local_locus.pl
    pdl_gunplot_local_gtf.pl
    select_top_hits_from_wide.pl
)) {
    my $path = File::Spec->catfile($canonical_dir, $file);
    die "Required gnuplot helper is missing: $path\n" unless -f $path;
}

my @active_files = (
    File::Spec->catfile($project_root, 'auto_prepare_and_run_diff_gwas_with_gunplot.pl'),
    File::Spec->catfile($Bin, 'generate_requested_top_hits_csv.pl'),
    File::Spec->catfile($project_root, 'install', 'check_pipeline_install.sh'),
);

for my $path (@active_files) {
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    local $/;
    my $text = <$fh>;
    close $fh or die "Cannot close $path: $!\n";

    die "Active executable path still names DiffGWASDeps/gunplot in $path\n"
        if $text =~ m{DiffGWASDeps[\\/]gunplot[\\/]};
    die "Active Perl path still uses the legacy gunplot directory in $path\n"
        if $text =~ /['\"]DiffGWASDeps['\"]\s*,\s*['\"]gunplot['\"]/;
}

my $wrapper = $active_files[0];
open my $wrapper_fh, '<', $wrapper or die "Cannot read $wrapper: $!\n";
local $/;
my $wrapper_text = <$wrapper_fh>;
close $wrapper_fh or die "Cannot close $wrapper: $!\n";

die "Canonical .gnuplot_gtf_cache default is missing from $wrapper\n"
    unless $wrapper_text =~ /\.gnuplot_gtf_cache/;
die "Legacy cache migration support is missing from $wrapper\n"
    unless $wrapper_text =~ /\.gunplot_gtf_cache/;

print "gnuplot directory layout ok\n";
