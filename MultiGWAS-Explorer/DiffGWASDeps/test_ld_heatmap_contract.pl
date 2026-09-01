#!/usr/bin/env perl
use strict;
use warnings;

use File::Spec;
use FindBin qw($Bin);

my %required = (
    'SNP_Local_Manhattan_With_GTF.sas' => [
        'LD_display_mode=none', 'LD_r2_values=',
        'ld_heatmap_var=&effective_ld_heatmap_var',
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

print "SAS/gnuplot optional LD heatmap contract: PASS\n";
