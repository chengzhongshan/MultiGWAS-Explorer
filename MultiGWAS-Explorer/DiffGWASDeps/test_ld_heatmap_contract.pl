#!/usr/bin/env perl
use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);

my %required = (
    'SNP_Local_Manhattan_With_GTF.sas' => [
        'LD_display_mode=none', 'LD_r2_values=',
        'ld_heatmap_var=&effective_ld_heatmap_var',
        'effective_signed_r2=1', '_SIGNED_R2_',
        'effective_heatmap_min_neg_val=-1',
        'effective_heatmap_max_pos_val=1',
    ],
    'Multgscatter_with_gene_exons.sas' => [
        'ld_heatmap_var=', 'ld_heatmap_var=&ld_heatmap_var',
    ],
    'map_grp_assoc2gene4covidsexgwas.sas' => [
        'ld_heatmap_var=', 'ld_heatmap_var=&ld_heatmap_var',
    ],
    'Lattice_gscatter_over_bed_track.sas' => [
        'ld_heatmap_var=', 'rangeattrmap name="ldheatmap"',
        'continuouslegend "ldsc"', 'heatmap_legend_title=',
    ],
    'Manhattan4DiffGWASs_png.sas' => [
        'LD_display_mode=none', '_LD_R2_=', 'LD_heatmap_legend_title=',
    ],
    'run_sas_oda_local_top_hits_with_gtf.sas' => [
        '__GTF_LD_DISPLAY_MODE__', '__GTF_LD_R2_VALUES__',
        '__GTF_LD_HEATMAP_COLORS__', '__GTF_LD_HEATMAP_LEGEND_TITLE__',
    ],
    'run_sas_oda_local_top_hits_manhattan.sas' => [
        '__GTF_LD_DISPLAY_MODE__', '__GTF_LD_R2_VALUES__',
        '__GTF_LD_HEATMAP_COLORS__', '__GTF_LD_HEATMAP_LEGEND_TITLE__',
    ],
);

for my $name (sort keys %required) {
    my $path = File::Spec->catfile($Bin, $name);
    open my $fh, '<:raw', $path or die "Cannot read $path: $!\n";
    my $text = do { local $/; <$fh> };
    close $fh;
    for my $token (@{ $required{$name} }) {
        die "$name is missing LD heatmap contract token: $token\n"
            unless index($text, $token) >= 0;
    }
}

for my $shell (qw(
    run_sas_oda_local_top_hits_with_gtf_download_html.sh
    run_sas_oda_local_top_hits_manhattan_download_png.sh
    run_sas_oda_single_snp_with_gtf_download_html.sh
)) {
    my $path = File::Spec->catfile($Bin, $shell);
    open my $fh, '<:raw', $path or die "Cannot read $path: $!\n";
    my $text = do { local $/; <$fh> };
    close $fh;
    for my $name (qw(GTF_LD_DISPLAY_MODE GTF_LD_R2_VALUES GTF_LD_HEATMAP_COLORS GTF_LD_HEATMAP_LEGEND_TITLE)) {
        die "$shell does not render $name\n"
            unless $text =~ /--replace\s+"$name=/;
    }
}

my $single_runner = File::Spec->catfile($Bin, 'run_sas_oda_single_snp_with_gtf_download_html.sh');
open my $single_fh, '<:raw', $single_runner or die "Cannot read $single_runner: $!\n";
my $single_text = do { local $/; <$single_fh> };
close $single_fh;
for my $token (
    'GTF_LD_R2_CACHE',
    'augment_gwas_with_ld_r2.pl',
    '--extra-numeric-cols LD_R2',
    'GTF_LD_REFERENCE_SNP',
) {
    die "Single-SNP SAS runner is missing numeric LD-cache support: $token\n"
        unless index($single_text, $token) >= 0;
}

my $root = File::Spec->catdir($Bin, File::Spec->updir());
my %reference_contract = (
    'auto_prepare_and_run_diff_gwas.pl' => [
        'local-ld-reference-snp|ld-reference-snp=s',
        'query_snp   => $local_ld_reference_snp',
        '_r2_${threshold_tag}_w${window_tag}.plink2_1kg_phase3.tsv',
        '1000G Phase 3 / PLINK2',
    ],
    'auto_prepare_and_run_diff_gwas_with_gunplot.pl' => [
        'ld-reference-snp=s',
        'query_snps => [$ld_reference_snp]',
        "push \@cmd, ('--ld-reference-snp', \$ld_reference_snp)",
    ],
    'server.pl' => [
        'local_ld_reference_snp',
        "('--local-ld-reference-snp', \$local_ld_reference_snp)",
        'ld_reference_snp',
        "('--ld-reference-snp', \$ld_reference_snp)",
    ],
);
for my $name (sort keys %reference_contract) {
    my $path = File::Spec->catfile($root, $name);
    open my $fh, '<:raw', $path or die "Cannot read $path: $!\n";
    my $text = do { local $/; <$fh> };
    close $fh;
    for my $token (@{ $reference_contract{$name} }) {
        die "$name is missing single-reference LD contract token: $token\n"
            unless index($text, $token) >= 0;
    }
}

my $sas_auto_path = File::Spec->catfile($root, 'auto_prepare_and_run_diff_gwas.pl');
open my $sas_auto_fh, '<:raw', $sas_auto_path or die "Cannot read $sas_auto_path: $!\n";
my $sas_auto_text = do { local $/; <$sas_auto_fh> };
close $sas_auto_fh;
for my $unsafe_label (
    'sign(Z); 1000 Genomes',
    '$local_ld_population_label; 1000G',
) {
    die "SAS label defaults contain an unquoted statement-ending semicolon: $unsafe_label\n"
        if index($sas_auto_text, $unsafe_label) >= 0;
}
for my $forbidden (
    'query_snps => join(\',\', @configured_target_snps)',
    'local_ld_cache_is_multi',
) {
    die "SAS local-GTF orchestration must retain one explicit LD reference per plotted locus; found: $forbidden\n"
        if index($sas_auto_text, $forbidden) >= 0;
}

print "SAS/gnuplot optional LD heatmap contract: PASS\n";
