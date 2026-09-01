#!/usr/bin/env perl
BEGIN {
    require File::Basename;
    require File::Spec;
    require Config;
    require Cwd;
    require lib;
    my $current_arch = lc($Config::Config{archname} || '');
    my $current_os = lc($^O || '');
    my $is_arch_dir = sub {
        my ($dir) = @_;
        return 0 unless -d $dir;
        my $name = File::Basename::basename($dir);
        my $looks_arch_specific = ($name =~ /(?:-thread-multi|linux|gnu|darwin|MSWin32|cygwin|^x86_64|^aarch64|^arm64|^i[3-6]86)/i)
            || -d File::Spec->catdir($dir, 'auto');
        return 0 unless $looks_arch_specific;
        return 1 if $current_arch && (lc($name) eq $current_arch || index($current_arch, lc($name)) >= 0 || index(lc($name), $current_arch) >= 0);
        return 1 if $current_os eq 'cygwin' && $name =~ /cygwin/i;
        return 1 if $current_os =~ /linux/  && $name =~ /(?:linux|gnu)/i;
        return 1 if $current_os =~ /darwin/ && $name =~ /darwin/i;
        return 1 if $current_os =~ /mswin32/ && $name =~ /MSWin32/i;
        return 0;
    };
    my $script_dir = Cwd::abs_path(File::Basename::dirname(__FILE__)) || File::Basename::dirname(__FILE__);
    my $platform_tag = lc($^O || '');
    $platform_tag =~ s/[^a-z0-9]+/_/g;
    for my $root ($script_dir, File::Spec->catdir($script_dir, File::Spec->updir())) {
        my @base_candidates;
        if (defined $ENV{PIPELINE_PERL_LOCAL_DIR} && length $ENV{PIPELINE_PERL_LOCAL_DIR}) {
            push @base_candidates, File::Spec->catdir($ENV{PIPELINE_PERL_LOCAL_DIR}, 'lib', 'perl5');
        }
        push @base_candidates, File::Spec->catdir($root, 'local', "perl5-$platform_tag", 'lib', 'perl5');
        push @base_candidates, File::Spec->catdir($root, 'local', 'perl5', 'lib', 'perl5');
        my %seen_base;
        for my $base (@base_candidates) {
            next unless -d $base;
            next if $seen_base{$base}++;
            lib->import($base);
            for my $arch (glob(File::Spec->catdir($base, '*'))) {
                next unless $is_arch_dir->($arch);
                lib->import($arch);
            }
        }
    }
}
use strict;
use warnings;

use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);
use JSON::PP qw(decode_json);
use GD;
use Text::ParseWords qw(parse_line);
use File::Spec;
use File::Basename qw(basename dirname);
use File::Copy qw(copy);
use File::Path qw(make_path);
use POSIX qw(strftime ceil);
use Time::HiRes qw(time);
use IO::Uncompress::Gunzip qw($GunzipError);

sub usage {
    return <<"USAGE";
Usage:
  perl auto_prepare_and_run_diff_gwas_with_gunplot.pl --spec spec.json [options]

Options:
  --plots LIST                  Comma list: manhattan,local_manhattan,local_gtf,forest
  --step NAME                   Plot step(s): plot_manhattan, plot_local_manhattan, plot_local_gtf, plot_forest
  --force                       Force upstream preprocessing refresh.
  --display-gwas LIST           Comma-separated displayed GWAS tracks, shared with the
                                SAS ODA runner config. Use pair prefixes such as
                                ALL,EUR,ASN for differential tracks and GWAS labels
                                such as ALL_FEMALE or EUR_MALE for single-GWAS tracks.
  --remove-X-chr / --no-remove-X-chr
                                Remove chromosome X from final gunplot figures.
                                Default: enabled.
  --target-snps A,B,C           Override target SNP list.
  --target-snp-genes MAP        Optional SNP:GENE overrides, comma-separated. Example: rs17425819:JAK2,rs2564978:CR1
  --ld-snps A,B,C               Optional LD-linked SNPs to overlay with asterisks.
  --[no-]highlight-high-ld-snps Resolve and draw high-LD markers. Default: off.
                                Supplying --ld-snps also enables highlighting.
  --ld-marker-symbol NAME       star, plus, cross, circle, square, triangle, or diamond.
                                Default: star when highlighting is enabled.
  --ld-marker-color COLOR       Named or #RRGGBB marker color. Default: black.
  --ld-display-mode MODE        none, markers, heatmap, or both. Default: none.
                                heatmap colors LD proxies by r2 and adds a separate
                                inset legend; it does not replace the Z-score scale.
  --ld-r2-values MAP            Optional SNP:r2 pairs, comma-separated, for explicit
                                --ld-snps (example rs1:0.92,rs2:0.81).
  --ld-heatmap-colors LIST      Low-to-high #RRGGBB colors for the LD inset.
                                Default: #f7fbff,#6baed6,#54278f.
  --ld-cache FILE               Normalized local HaploReg TSV/SQLite LD cache.
  --ld-population POP           HaploReg population for high-LD stars (default EUR).
  --ld-r2-threshold N           High-LD star threshold (default 0.8).
  --[no-]ld-web-fallback        Query HaploReg only when the local cache misses
                                a query SNP (default enabled); successful queries
                                are saved in cache/haploreg_ld for reuse.
  --ld-audit-file FILE          Optional LD audit TSV; PRUNED_LD candidates linked
                                to a query lead SNP are marked automatically.
  --get-common-associations [VALUE]
                                Forward common-association mode to upstream preprocessing.
  --common-association-top-hit-threshold VALUE
                                Optional starting threshold for common-association mode.
  --reference-build BUILD      Override genome build for local GTF annotations:
                                hg19, hg38, or t2t. Default: auto-detect, then hg38 fallback.
  --local-gtf-window-bp N       Override local GTF window.
  --local-max-hits-per-fig N    Override local Manhattan batch size.
  --local-manhattan-columns N   Override the number of loci columns per combined gunplot local-Manhattan figure.
  --local-manhattan-annotation MODE  Under-column annotation for combined local Manhattan: labels, gtf, auto, none
                                     Default: gtf

Recommended interactive pattern:
  open this repository in VS Code, use the integrated shell for the current
  platform, and keep Codex in the same workspace if you are driving the local
  gunplot pipeline through the MCP server.
USAGE
}

my $spec_file = '';
my $plots = 'manhattan,local_manhattan,local_gtf';
my @step_args;
my $force = 0;
my $display_gwas_override = '';
my $target_snps_override = '';
my $target_snp_genes_override = '';
my $ld_snps_override = '';
my $ld_audit_file_override = '';
my $ld_cache_override = '';
my $ld_population_override = 'EUR';
my $ld_r2_threshold_override = 0.8;
my $ld_web_fallback = 1;
my $highlight_high_ld_snps = 0;
my $ld_marker_symbol = 'star';
my $ld_marker_color = 'black';
my $ld_display_mode = 'none';
my $ld_r2_values_override = '';
my $ld_heatmap_colors = '#f7fbff,#6baed6,#54278f';
my $get_common_associations = '';
my $common_association_top_hit_threshold = '';
my $reference_build_override = '';
my $local_gtf_window_bp_override = '';
my $local_max_hits_per_fig_override = 0;
my $local_manhattan_columns_override = 0;
my $local_manhattan_annotation_override = '';
my $remove_x_chr = 1;
my $help = 0;

GetOptions(
    'spec=s'                  => \$spec_file,
    'plots=s'                 => \$plots,
    'step=s@'                 => \@step_args,
    'force!'                  => \$force,
    'display-gwas|display-tracks=s' => \$display_gwas_override,
    'remove-x-chr!'            => \$remove_x_chr,
    'target-snps=s'           => \$target_snps_override,
    'target-snp-genes=s'      => \$target_snp_genes_override,
    'ld-snps=s'               => \$ld_snps_override,
    'ld-audit-file=s'         => \$ld_audit_file_override,
    'ld-cache=s'              => \$ld_cache_override,
    'ld-population=s'         => \$ld_population_override,
    'ld-r2-threshold=f'       => \$ld_r2_threshold_override,
    'ld-web-fallback!'        => \$ld_web_fallback,
    'highlight-high-ld-snps!' => \$highlight_high_ld_snps,
    'ld-marker-symbol=s'      => \$ld_marker_symbol,
    'ld-marker-color=s'       => \$ld_marker_color,
    'ld-display-mode=s'       => \$ld_display_mode,
    'ld-r2-values=s'          => \$ld_r2_values_override,
    'ld-heatmap-colors=s'     => \$ld_heatmap_colors,
    'get-common-associations:s' => \$get_common_associations,
    'common-association-top-hit-threshold=s' => \$common_association_top_hit_threshold,
    'reference-build=s'       => \$reference_build_override,
    'local-gtf-window-bp=s'   => \$local_gtf_window_bp_override,
    'local-max-hits-per-fig=i'=> \$local_max_hits_per_fig_override,
    'local-manhattan-columns=i'=> \$local_manhattan_columns_override,
    'local-manhattan-annotation=s'=> \$local_manhattan_annotation_override,
    'help!'                   => \$help,
) or die usage();

if ($help || !$spec_file) {
    print usage();
    exit($help ? 0 : 1);
}
$ld_population_override = uc($ld_population_override || 'EUR');
die "--ld-population must be AFR, AMR, ASN, or EUR\n"
    unless $ld_population_override =~ /^(?:AFR|AMR|ASN|EUR)$/;
die "--ld-r2-threshold must be between 0 and 1\n"
    unless $ld_r2_threshold_override > 0 && $ld_r2_threshold_override <= 1;
$ld_display_mode = lc(trim($ld_display_mode || 'none'));
die "--ld-display-mode must be none, markers, heatmap, or both\n"
    unless $ld_display_mode =~ /^(?:none|markers|heatmap|both)$/;
$ld_display_mode = 'markers'
    if $ld_display_mode eq 'none' && ($highlight_high_ld_snps || length(trim($ld_snps_override)));
$highlight_high_ld_snps = 1 if $ld_display_mode ne 'none';
$ld_marker_symbol = lc(trim($ld_marker_symbol || 'star'));
die "--ld-marker-symbol must be star, plus, cross, circle, square, triangle, or diamond\n"
    unless $ld_marker_symbol =~ /^(?:star|plus|cross|circle|square|triangle|diamond)$/;
die "--ld-marker-color must be a named color or #RRGGBB\n"
    unless $ld_marker_color =~ /^(?:#[0-9A-Fa-f]{6}|[A-Za-z][A-Za-z0-9_-]*)$/;
my @ld_heatmap_colors = split /,/, $ld_heatmap_colors;
die "--ld-heatmap-colors requires at least two comma-separated #RRGGBB colors\n"
    unless @ld_heatmap_colors >= 2
        && !grep { trim($_) !~ /^#[0-9A-Fa-f]{6}$/ } @ld_heatmap_colors;
$ld_heatmap_colors = join(',', map { lc(trim($_)) } @ld_heatmap_colors);

die "Spec file not found: $spec_file\n" unless -f $spec_file;

my %requested = normalize_requested_plots($plots, \@step_args);
die "No gunplot plot steps were requested.\n" unless grep { $requested{$_} } qw(plot_manhattan plot_local_manhattan plot_local_gtf plot_forest);

my $spec = load_json($spec_file);
my $target_snp_gene_overrides = parse_target_snp_gene_overrides(
    spec_value => $spec->{target_snp_genes},
    cli_value  => $target_snp_genes_override,
);
my $artifact_stem = $spec->{artifact_stem} || $spec->{project_tag} || 'diff_gwas';
my $workdir_local = resolve_portable_workdir($spec->{workdir} || $Bin, $Bin);
my $configs_dir_local = File::Spec->catdir($workdir_local, 'configs');
my $output_dir_local = localize_path($spec->{output_dir} || $spec->{input_dir} || $Bin);
my $shared_gtf_cache_local = localize_path(
    $spec->{gtf_cache_dir} || File::Spec->catdir($workdir_local, 'cache', 'gtf')
);

make_path($configs_dir_local) unless -d $configs_dir_local;
make_path($output_dir_local) unless -d $output_dir_local;
make_path($shared_gtf_cache_local) unless -d $shared_gtf_cache_local;

my $runner_config_local = File::Spec->catfile($configs_dir_local, "auto_${artifact_stem}_runner.json");
my $preset_config_local = File::Spec->catfile($configs_dir_local, "auto_${artifact_stem}_preset.json");

my $reused_existing_runner = 0;
my $has_runner_override = 0;
for my $override_value (
    $display_gwas_override,
    $target_snps_override,
    $target_snp_genes_override,
    $get_common_associations,
    $common_association_top_hit_threshold,
    $reference_build_override,
    $local_gtf_window_bp_override,
) {
    if (defined $override_value && length $override_value) {
        $has_runner_override = 1;
        last;
    }
}
$has_runner_override = 1 if $local_max_hits_per_fig_override;

if (!$force && !$has_runner_override && -f $runner_config_local) {
    my $existing_runner = load_json($runner_config_local);
    my $existing_data = localize_path($existing_runner->{DATA_GZ} || '');
    if ($existing_data && -s $existing_data) {
        my ($ok, $why_not) = runner_wide_cache_is_reusable($existing_runner);
        if ($ok) {
        $reused_existing_runner = 1;
        print "[skip] reusing existing runner config and wide subset for gunplot pipeline\n";
        }
        else {
            print "[info] refreshing upstream preprocessing because the cached wide subset is not reusable: $why_not\n";
        }
    }
}

if (!$reused_existing_runner) {
    run_upstream_preprocessing(
        spec_file                       => $spec_file,
        force                           => $force,
        display_gwas_override           => $display_gwas_override,
        target_snps_override            => $target_snps_override,
        target_snp_genes_override       => $target_snp_genes_override,
        get_common_associations         => $get_common_associations,
        common_association_top_hit_threshold => $common_association_top_hit_threshold,
        reference_build_override        => $reference_build_override,
        local_gtf_window_bp_override    => $local_gtf_window_bp_override,
        local_max_hits_per_fig_override => $local_max_hits_per_fig_override,
    );
}

die "Runner config was not generated: $runner_config_local\n" unless -f $runner_config_local;
my $runner = load_json($runner_config_local);

my $gnuplot = find_gnuplot_exe();
print "Using gnuplot executable: $gnuplot\n";

my $wide_data_local = localize_path($runner->{DATA_GZ});
die "Wide input file not found: $wide_data_local\n" unless -s $wide_data_local;

my @manhattan_pcols = (
    $runner->{MANHATTAN_P_VAR},
    @{ ref($runner->{MANHATTAN_OTHER_P_VARS}) eq 'ARRAY' ? $runner->{MANHATTAN_OTHER_P_VARS} : [] },
);

my @manhattan_labels = split /\|/, ($runner->{MANHATTAN_GWAS_LABEL_NAMES} || join('|', @manhattan_pcols));

my @gtf_pcols = grep { length } split /\s+/, ($runner->{GTF_ASSOC_PVARS} || '');
my @gtf_zcols = grep { length } split /\s+/, ($runner->{GTF_ZSCORE_VARS} || '');
my @gtf_labels = grep { length } split /\s+/, ($runner->{GTF_LABELS} || '');
@gtf_labels = @gtf_pcols unless @gtf_labels == @gtf_pcols;

my $indexed_source_long_local = '';
if ($requested{plot_local_manhattan} || $requested{plot_local_gtf}) {
    $indexed_source_long_local = ensure_tabix_ready_long_source(
        source_long => localize_path($runner->{SOURCE_LONG_GZ} || ''),
        output_dir  => $output_dir_local,
        gnuplot     => $gnuplot,
        force       => $force,
    );
}

my %outputs;
my @new_manhattan_labels=@manhattan_labels;
   @new_manhattan_labels=map{$_=~s/_/ /g;$_}@new_manhattan_labels;
if ($requested{plot_manhattan}) {
    my $step_started = time();
    %outputs = (
        %outputs,
        plot_manhattan(
            output_dir => $output_dir_local,
            runner     => $runner,
            wide_data  => $wide_data_local,
            gnuplot    => $gnuplot,
            pcols      => \@manhattan_pcols,
            labels     => \@new_manhattan_labels,
            remove_x_chr => $remove_x_chr,
            force      => $force,
        ),
    );
    print "[done] plot_manhattan finished in " . format_elapsed_seconds(time() - $step_started) . "\n";
}

my @hits;
if ($requested{plot_local_manhattan} || $requested{plot_local_gtf}) {
    my $effective_target_snps = $target_snps_override || ($runner->{TARGET_SNP_LIST} || '');
    @hits = collect_top_hits(
        spec_file    => $spec_file,
        runner_config_path => $runner_config_local,
        output_dir   => $output_dir_local,
        runner       => $runner,
        wide_data    => $wide_data_local,
        target_snps  => $effective_target_snps,
        remove_x_chr => $remove_x_chr,
        force        => $force,
    );
    enrich_hits_from_local_top_hits_csv(
        output_dir => $output_dir_local,
        runner     => $runner,
        hits       => \@hits,
    );
    apply_target_snp_gene_overrides(
        hits      => \@hits,
        overrides => $target_snp_gene_overrides,
    );
    if ($remove_x_chr && !$effective_target_snps) {
        @hits = grep { !is_chr_x($_->{CHR}) } @hits;
    }
    die "No local top hits were selected.\n" unless @hits;
}

if ($requested{plot_local_manhattan}) {
    my $step_started = time();
    %outputs = (
        %outputs,
        plot_local_series(
            kind         => 'local_manhattan',
            output_dir   => $output_dir_local,
            runner       => $runner,
            wide_data    => $wide_data_local,
            gnuplot      => $gnuplot,
            hits         => \@hits,
            pcols        => \@manhattan_pcols,
            labels       => \@manhattan_labels,
            window_bp    => ($runner->{LOCAL_WINDOW_BP} || '1e7'),
            batch_size   => ($local_max_hits_per_fig_override || ($runner->{LOCAL_MAX_HITS_PER_FIG} || 15)),
            panel_columns=> ($local_manhattan_columns_override || $runner->{LOCAL_MANHATTAN_COLUMNS} || $spec->{local_manhattan_columns} || 0),
            annotation_mode => ($local_manhattan_annotation_override
                || $runner->{LOCAL_MANHATTAN_COLUMN_ANNOTATION}
                || $spec->{local_manhattan_annotation}
                || 'gtf'),
            output_base  => gunplotize_name($runner->{LOCAL_OUTPUT_PREFIX} || 'local_top_hits_manhattan'),
            html_title   => gunplot_title($runner->{LOCAL_HTML_TITLE} || 'Local top hits Manhattan Plot'),
            height       => ($runner->{LOCAL_MANHATTAN_FIG_HEIGHT} || 1200),
            with_gtf     => 0,
            source_long  => ($indexed_source_long_local || localize_path($runner->{SOURCE_LONG_GZ} || '')),
            preset_config=> $preset_config_local,
            gtf_cache_dir => $shared_gtf_cache_local,
            top_csv_name => gunplotize_name($runner->{LOCAL_TOP_HITS_CSV_BASENAME} || 'gunplot_top_hits.csv'),
            ld_snps      => $ld_snps_override,
            ld_audit_file=> $ld_audit_file_override,
            ld_cache     => $ld_cache_override,
            ld_population=> $ld_population_override,
            ld_r2_threshold => $ld_r2_threshold_override,
            ld_web_fallback => $ld_web_fallback,
            highlight_high_ld_snps => $highlight_high_ld_snps,
            ld_marker_symbol => $ld_marker_symbol,
            ld_marker_color  => $ld_marker_color,
            ld_display_mode  => $ld_display_mode,
            ld_r2_values     => $ld_r2_values_override,
            ld_heatmap_colors=> $ld_heatmap_colors,
            force       => $force,
        ),
    );
    print "[done] plot_local_manhattan finished in " . format_elapsed_seconds(time() - $step_started) . "\n";
}

if ($requested{plot_local_gtf}) {
    my $gtf_window = $local_gtf_window_bp_override || ($runner->{LOCAL_GTF_WINDOW_BP} || $runner->{LOCAL_WINDOW_BP} || '1e7');
    my $step_started = time();
    %outputs = (
        %outputs,
        plot_local_series(
            kind         => 'local_gtf',
            output_dir   => $output_dir_local,
            runner       => $runner,
            wide_data    => $wide_data_local,
            gnuplot      => $gnuplot,
            hits         => \@hits,
            pcols        => \@gtf_pcols,
            zcols        => \@gtf_zcols,
            labels       => \@gtf_labels,
            window_bp    => $gtf_window,
            batch_size   => ($runner->{LOCAL_GTF_MAX_HITS_PER_FIG} || 1),
            output_base  => gunplotize_name($runner->{OUTPUT_HTML_BASENAME} || 'local_top_hits_with_gtf.html'),
            html_title   => gunplot_title($runner->{LOCAL_HTML_TITLE} || 'Local top hits Manhattan and GTF Plot'),
            height       => compute_gtf_height(scalar(@gtf_pcols)),
            with_gtf     => 1,
            source_long  => ($indexed_source_long_local || localize_path($runner->{SOURCE_LONG_GZ} || '')),
            preset_config=> $preset_config_local,
            gtf_cache_dir => $shared_gtf_cache_local,
            top_csv_name => gunplotize_name($runner->{LOCAL_TOP_HITS_CSV_BASENAME} || 'gunplot_top_hits.csv'),
            ld_snps      => $ld_snps_override,
            ld_audit_file=> $ld_audit_file_override,
            ld_cache     => $ld_cache_override,
            ld_population=> $ld_population_override,
            ld_r2_threshold => $ld_r2_threshold_override,
            ld_web_fallback => $ld_web_fallback,
            highlight_high_ld_snps => $highlight_high_ld_snps,
            ld_marker_symbol => $ld_marker_symbol,
            ld_marker_color  => $ld_marker_color,
            ld_display_mode  => $ld_display_mode,
            ld_r2_values     => $ld_r2_values_override,
            ld_heatmap_colors=> $ld_heatmap_colors,
            force       => $force,
        ),
    );
    print "[done] plot_local_gtf finished in " . format_elapsed_seconds(time() - $step_started) . "\n";
}

if ($requested{plot_forest}) {
    my $step_started = time();
    %outputs = (
        %outputs,
        plot_forest(
            spec_file    => $spec_file,
            output_dir   => $output_dir_local,
            runner       => $runner,
            runner_config=> $runner_config_local,
            wide_data    => $wide_data_local,
            gnuplot      => $gnuplot,
            target_snps  => ($target_snps_override || ($runner->{TARGET_SNP_LIST} || '')),
            target_snp_genes => ($target_snp_genes_override || ($runner->{TARGET_SNP_GENES} || '')),
            remove_x_chr => $remove_x_chr,
            force        => $force,
        ),
    );
    print "[done] plot_forest finished in " . format_elapsed_seconds(time() - $step_started) . "\n";
}

print "\nGenerated gunplot outputs:\n";
for my $key (sort keys %outputs) {
    print "$key\t$outputs{$key}\n";
}

exit 0;

sub run_upstream_preprocessing {
    my (%args) = @_;
    my @cmd = (
        $^X,
        File::Spec->catfile($Bin, 'auto_prepare_and_run_diff_gwas.pl'),
        '--spec', $args{spec_file},
        '--step', 'extract_wide_subset',
    );
    push @cmd, '--force' if $args{force};
    if ($args{display_gwas_override}) {
        push @cmd, '--display-gwas', $args{display_gwas_override};
    }
    if ($args{target_snps_override}) {
        push @cmd, '--target-snps', $args{target_snps_override};
    }
    if ($args{target_snp_genes_override}) {
        push @cmd, '--target-snp-genes', $args{target_snp_genes_override};
    }
    if (defined $args{get_common_associations} && length $args{get_common_associations}) {
        if ($args{get_common_associations} =~ /^(?:1|true|yes|y)$/i) {
            push @cmd, '--get-common-associations';
        }
        else {
            push @cmd, "--get-common-associations=$args{get_common_associations}";
        }
    }
    if ($args{common_association_top_hit_threshold}) {
        push @cmd, '--common-association-top-hit-threshold', $args{common_association_top_hit_threshold};
    }
    if ($args{reference_build_override}) {
        push @cmd, '--reference-build', $args{reference_build_override};
    }
    if ($args{local_gtf_window_bp_override}) {
        push @cmd, '--local-gtf-window-bp', $args{local_gtf_window_bp_override};
    }
    if ($args{local_max_hits_per_fig_override}) {
        push @cmd, '--local-max-hits-per-fig', $args{local_max_hits_per_fig_override};
    }
    run_cmd(\@cmd, 'upstream preprocessing');
}

sub plot_manhattan {
    my (%args) = @_;
    my $output_prefix = gunplotize_name($args{runner}{OUTPUT_PREFIX} || 'gunplot_manhattan');
    my $out_prefix_path = File::Spec->catfile($args{output_dir}, $output_prefix);
    my $png_path = $out_prefix_path . '.png';
    my $html = File::Spec->catfile($args{output_dir}, $output_prefix . '.html');
    my $min_logp = defined $args{runner}{GUNPLOT_MANHATTAN_MIN_LOGP}
        ? $args{runner}{GUNPLOT_MANHATTAN_MIN_LOGP}
        : 0.5;
    my $keep_all_logp = defined $args{runner}{GUNPLOT_MANHATTAN_KEEP_ALL_LOGP}
        ? $args{runner}{GUNPLOT_MANHATTAN_KEEP_ALL_LOGP}
        : 2;
    my $thin_mod = defined $args{runner}{GUNPLOT_MANHATTAN_THIN_MOD}
        ? $args{runner}{GUNPLOT_MANHATTAN_THIN_MOD}
        : 10;
    if (!(!$args{force} && -s $png_path)) {
        my @cmd = (
            $^X,
            File::Spec->catfile($Bin, 'DiffGWASDeps', 'gnuplot', 'pdl_gunplot_manhattan.pl'),
            '--data', $args{wide_data},
            '--out-prefix', $out_prefix_path,
            '--pcols', join(',', @{ $args{pcols} || [] }),
            '--labels', join('|', @{ $args{labels} || [] }),
            '--width', ($args{runner}{MANHATTAN_FIG_WIDTH} || 1800),
            '--height', ($args{runner}{MANHATTAN_FIG_HEIGHT} || 1250),
            '--sig', ($args{runner}{TOP_HIT_SIGNAL_THRSHD} || '1e-6'),
            '--min-logp', $min_logp,
            '--keep-all-logp', $keep_all_logp,
            '--thin-mod', $thin_mod,
            '--gnuplot', $args{gnuplot},
        );
        push @cmd, '--remove-x-chr' if $args{remove_x_chr};
        run_cmd(\@cmd, 'gunplot genomewide Manhattan');
    }
    else {
        print "[skip] reusing existing gunplot genomewide Manhattan PNG $png_path\n";
    }

    write_single_image_html(
        path      => $html,
        title     => gunplot_title($args{runner}{HTML_TITLE} || 'Gunplot Manhattan Plot'),
        image_rel => basename($png_path),
    );
    return (
        gunplot_manhattan_png  => $png_path,
        gunplot_manhattan_html => $html,
    );
}

sub plot_forest {
    my (%args) = @_;
    my $runner = $args{runner} || {};
    my $forest_output_prefix = gunplotize_name($runner->{FOREST_OUTPUT_PREFIX} || 'gunplot_top_hits_forest');
    $forest_output_prefix =~ s/\.html$//i;
    my $out_prefix_path = File::Spec->catfile($args{output_dir}, $forest_output_prefix);
    my $html_name = gunplotize_name($runner->{FOREST_OUTPUT_HTML_BASENAME} || ($forest_output_prefix . '.html'));
    my $html_path = File::Spec->catfile($args{output_dir}, $html_name);
    my $manifest_name = gunplotize_name($runner->{FOREST_OUTPUT_MANIFEST_BASENAME} || ($forest_output_prefix . '.manifest.tsv'));
    my $manifest_path = File::Spec->catfile($args{output_dir}, $manifest_name);
    my $top_hits_csv_name = forest_csv_basename_for_targets(
        base_name    => ($runner->{FOREST_TOP_HITS_CSV_BASENAME} || ($forest_output_prefix . '_top_hits.csv')),
        target_snps  => $args{target_snps},
    );
    $top_hits_csv_name = gunplotize_name($top_hits_csv_name);
    my $top_hits_csv_path = File::Spec->catfile($args{output_dir}, $top_hits_csv_name);

    if ($args{force} || !-s $top_hits_csv_path) {
        my @cmd = (
            $^X,
            File::Spec->catfile($Bin, 'DiffGWASDeps', 'generate_requested_top_hits_csv.pl'),
            '--input', $args{wide_data},
            '--output', $top_hits_csv_path,
            '--runner-config', $args{runner_config},
            '--top-hit-mode', ($runner->{TOP_HIT_MODE} || 'differential'),
            '--top-hit-focus-pvar', ($runner->{TOP_HIT_FOCUS_PVAR} || $runner->{MANHATTAN_P_VAR} || 'ALL_STD_P'),
            '--top-hit-signal-thrshd', ($runner->{TOP_HIT_SIGNAL_THRSHD} || '1e-6'),
            '--top-hit-signal-thrshds', ($runner->{TOP_HIT_SIGNAL_THRSHDS} || $runner->{TOP_HIT_SIGNAL_THRSHD} || '1e-6'),
            '--top-hit-dist-bp', ($runner->{TOP_HIT_DIST_BP} || '1e6'),
            '--maf-threshold', runner_maf_threshold($runner),
            '--max-hits', (defined $runner->{TOP_HIT_MAX_LOCI} ? $runner->{TOP_HIT_MAX_LOCI} : 0),
        );
        push @cmd, ('--spec', $args{spec_file}) if defined $args{spec_file} && length $args{spec_file};
        push @cmd, ('--target-snps', $args{target_snps}) if defined $args{target_snps} && length $args{target_snps};
        push @cmd, ('--target-snp-genes', $args{target_snp_genes}) if defined $args{target_snp_genes} && length $args{target_snp_genes};
        push @cmd, ('--gnomad-freq-file', $runner->{TOP_HIT_GNOMAD_FREQ_FILE})
            if defined $runner->{TOP_HIT_GNOMAD_FREQ_FILE} && length $runner->{TOP_HIT_GNOMAD_FREQ_FILE};
        push @cmd, ('--gnomad-pop-map', $runner->{TOP_HIT_GNOMAD_POP_MAP})
            if defined $runner->{TOP_HIT_GNOMAD_POP_MAP} && length $runner->{TOP_HIT_GNOMAD_POP_MAP};
        push @cmd, '--remove-x-chr' if $args{remove_x_chr};
        run_cmd(\@cmd, 'gunplot forest top-hit CSV generation');
    }
    else {
        print "[skip] reusing existing gunplot forest top-hit CSV $top_hits_csv_path\n";
    }

    my @png_candidates = (
        $out_prefix_path . '_single_snp.png',
        map { $out_prefix_path . '_' . safe_name($_) . '.png' } split(/\|/, ($runner->{FOREST_TRACK_IDS} || '')),
    );
    my $need_render = $args{force} || !-s $manifest_path;
    if (!$need_render) {
        my @existing_pngs = grep { -s $_ } @png_candidates;
        $need_render = @existing_pngs ? 0 : 1;
    }

    if ($need_render) {
        my @cmd = (
            $^X,
            File::Spec->catfile($Bin, 'DiffGWASDeps', 'gnuplot', 'pdl_gunplot_forest.pl'),
            '--csv', $top_hits_csv_path,
            '--out-prefix', $out_prefix_path,
            '--track-ids', ($runner->{FOREST_TRACK_IDS} || ''),
            '--track-labels', ($runner->{FOREST_TRACK_LABELS} || ''),
            '--track-beta-vars', ($runner->{FOREST_TRACK_BETA_VARS} || ''),
            '--track-se-vars', ($runner->{FOREST_TRACK_SE_VARS} || ''),
            '--track-p-vars', ($runner->{FOREST_TRACK_P_VARS} || ''),
            '--title', gunplot_title($runner->{FOREST_HTML_TITLE} || 'Top-hit forest plots'),
            '--width', ($runner->{FOREST_FIG_WIDTH} || 1100),
            '--dotsize', ($runner->{FOREST_DOTSIZE} || 8),
            '--y-font-size', ($runner->{FOREST_Y_FONT_SIZE} || 12),
            '--min-axis', ($runner->{FOREST_MIN_AXIS} || 0.4),
            '--max-axis', ($runner->{FOREST_MAX_AXIS} || 1.6),
            '--xaxis-ticks', ($runner->{FOREST_XAXIS_VALUE_RANGE} || '0.4 to 1.6 by 0.2'),
            '--default-hit-class', ($runner->{FOREST_DEFAULT_HIT_CLASS} || 'DIFFERENTIAL'),
            '--gnuplot', $args{gnuplot},
        );
        push @cmd, ('--height', $runner->{FOREST_FIG_HEIGHT})
            if defined $runner->{FOREST_FIG_HEIGHT} && length $runner->{FOREST_FIG_HEIGHT};
        run_cmd(\@cmd, 'gunplot forest plot rendering');
    }
    else {
        print "[skip] reusing existing gunplot forest plot artifacts rooted at $out_prefix_path\n";
    }

    my @manifest_rows = read_forest_manifest_rows($manifest_path);
    die "Gunplot forest manifest has no panel rows: $manifest_path\n" unless @manifest_rows;
    my @images = map {
        +{
            snp   => ($_->{track_label} || $_->{track_id} || 'forest'),
            image => File::Spec->catfile($args{output_dir}, $_->{png_file}),
        }
    } grep { defined $_->{png_file} && length $_->{png_file} && -s File::Spec->catfile($args{output_dir}, $_->{png_file}) } @manifest_rows;
    die "Gunplot forest manifest did not resolve any downloaded PNG panels under $args{output_dir}\n" unless @images;

    my $combined_png = '';
    if (@images > 1) {
        $combined_png = $out_prefix_path . '_combined.png';
        if ($args{force} || !png_is_newer_than_dependencies($combined_png, map { $_->{image} } @images)) {
            compose_png_grid(
                output_png => $combined_png,
                images     => [ map { $_->{image} } @images ],
                columns    => scalar(@images),
            );
        }
    }

    if (@images == 1) {
        write_single_image_html(
            path      => $html_path,
            title     => gunplot_title($runner->{FOREST_HTML_TITLE} || 'Top-hit forest plot'),
            image_rel => basename($images[0]{image}),
        );
    }
    elsif ($combined_png && -s $combined_png) {
        write_single_image_html(
            path      => $html_path,
            title     => gunplot_title($runner->{FOREST_HTML_TITLE} || 'Top-hit forest plots'),
            image_rel => basename($combined_png),
        );
    }
    else {
        write_gallery_html(
            path   => $html_path,
            title  => gunplot_title($runner->{FOREST_HTML_TITLE} || 'Top-hit forest plots'),
            images => \@images,
            top_csv=> basename($top_hits_csv_path),
        );
    }

    my %out = (
        gunplot_forest_html => $html_path,
        gunplot_forest_manifest => $manifest_path,
        gunplot_forest_top_hits_csv => $top_hits_csv_path,
    );
    if (@images == 1) {
        $out{gunplot_forest_png} = $images[0]{image};
    }
    elsif ($combined_png && -s $combined_png) {
        $out{gunplot_forest_png} = $combined_png;
        $out{gunplot_forest_panels} = scalar(@images);
    }
    else {
        $out{gunplot_forest_panels} = scalar(@images);
    }
    return %out;
}

sub collect_top_hits {
    my (%args) = @_;
    my $runner = $args{runner};
    my $target_snps = $args{target_snps} || '';
    my $differential_thresholds = resolve_runner_differential_threshold_ladder($runner);
    my $out_tsv = File::Spec->catfile(
        $args{output_dir},
        gunplotize_name(($runner->{PROJECT_TAG} || 'diff_gwas') . '.gunplot_top_hits.tsv')
    );

    if (!$args{force} && -s $out_tsv && !$target_snps) {
        print "[skip] reusing existing gunplot top-hit table $out_tsv\n";
        return read_hits_tsv($out_tsv);
    }

    if ($target_snps) {
        my @hits = map { +{ SNP => $_ } } grep { length } map { trim($_) } split /,/, $target_snps;
        return @hits;
    }

    if (($runner->{TOP_HIT_MODE} || '') =~ /^common_association$/i) {
        my @cmd = (
            $^X,
            File::Spec->catfile($Bin, 'DiffGWASDeps', 'verify_common_association_loci.pl'),
            '--spec', $args{spec_file},
            '--input', $args{wide_data},
            '--output', $out_tsv,
            '--top-p-thresholds', ($runner->{TOP_HIT_SIGNAL_THRSHDS} || $runner->{TOP_HIT_SIGNAL_THRSHD} || '1e-6'),
            '--top-hit-dist-bp', ($runner->{TOP_HIT_DIST_BP} || '1e8'),
            '--max-loci', 15,
            '--maf-threshold', runner_maf_threshold($runner),
        );
        push @cmd, ('--runner-config', $args{runner_config_path})
            if defined $args{runner_config_path} && length $args{runner_config_path};
        push @cmd, ('--gnomad-freq-file', $runner->{TOP_HIT_GNOMAD_FREQ_FILE})
            if defined $runner->{TOP_HIT_GNOMAD_FREQ_FILE} && length $runner->{TOP_HIT_GNOMAD_FREQ_FILE};
        push @cmd, ('--gnomad-pop-map', $runner->{TOP_HIT_GNOMAD_POP_MAP})
            if defined $runner->{TOP_HIT_GNOMAD_POP_MAP} && length $runner->{TOP_HIT_GNOMAD_POP_MAP};
        push @cmd, '--remove-x-chr' if $args{remove_x_chr};
        run_cmd(\@cmd, 'common-association top-hit selection');
        return read_hits_tsv($out_tsv);
    }
    if (($runner->{TOP_HIT_MODE} || '') =~ /^(?:common_and_differential|union_common_and_differential)$/i) {
        my @diff = collect_top_hits_for_mode(
            mode => 'differential',
            spec_file => $args{spec_file},
            output_dir => $args{output_dir},
            runner => $runner,
            runner_config_path => $args{runner_config_path},
            wide_data => $args{wide_data},
            remove_x_chr => $args{remove_x_chr},
        );
        my @common = collect_top_hits_for_mode(
            mode => 'common_association',
            spec_file => $args{spec_file},
            output_dir => $args{output_dir},
            runner => $runner,
            runner_config_path => $args{runner_config_path},
            wide_data => $args{wide_data},
            remove_x_chr => $args{remove_x_chr},
        );
        return merge_hit_lists_for_gunplot(
            max_hits => 15,
            lists => [ \@diff, \@common ],
        );
    }

    return collect_top_hits_for_mode(
        mode => 'differential',
        spec_file => $args{spec_file},
        output_dir => $args{output_dir},
        runner => $runner,
        runner_config_path => $args{runner_config_path},
        wide_data => $args{wide_data},
        remove_x_chr => $args{remove_x_chr},
        output_tsv => $out_tsv,
    );
}

sub runner_maf_threshold {
    my ($runner) = @_;
    return 0.01 unless ref($runner) eq 'HASH' && exists $runner->{TOP_HIT_MAF_THRESHOLD};
    my $value = $runner->{TOP_HIT_MAF_THRESHOLD};
    return 0 if defined $value && $value =~ /^(?:0|0\.0+)?$/;
    return $value if defined $value && $value ne '';
    return 0.01;
}

sub collect_top_hits_for_mode {
    my (%args) = @_;
    my $runner = $args{runner} || {};
    my $mode = $args{mode} || 'differential';
    my $out_tsv = $args{output_tsv}
        || File::Spec->catfile($args{output_dir}, "gunplot_top_hits_${mode}.tsv");
    if ($mode =~ /^common_association$/i) {
        my @cmd = (
            $^X,
            File::Spec->catfile($Bin, 'DiffGWASDeps', 'verify_common_association_loci.pl'),
            '--spec', $args{spec_file},
            '--input', $args{wide_data},
            '--output', $out_tsv,
            '--top-p-thresholds', ($runner->{TOP_HIT_SIGNAL_THRSHDS} || $runner->{TOP_HIT_SIGNAL_THRSHD} || '1e-6'),
            '--top-hit-dist-bp', ($runner->{TOP_HIT_DIST_BP} || '1e8'),
            '--max-loci', 15,
            '--maf-threshold', runner_maf_threshold($runner),
        );
        push @cmd, ('--runner-config', $args{runner_config_path})
            if defined $args{runner_config_path} && length $args{runner_config_path};
        push @cmd, ('--gnomad-freq-file', $runner->{TOP_HIT_GNOMAD_FREQ_FILE})
            if defined $runner->{TOP_HIT_GNOMAD_FREQ_FILE} && length $runner->{TOP_HIT_GNOMAD_FREQ_FILE};
        push @cmd, ('--gnomad-pop-map', $runner->{TOP_HIT_GNOMAD_POP_MAP})
            if defined $runner->{TOP_HIT_GNOMAD_POP_MAP} && length $runner->{TOP_HIT_GNOMAD_POP_MAP};
        push @cmd, '--remove-x-chr' if $args{remove_x_chr};
        run_cmd(\@cmd, 'common-association top-hit selection');
        return read_hits_tsv($out_tsv);
    }
    my $differential_thresholds = resolve_runner_differential_threshold_ladder($runner);
    my @cmd = (
        $^X,
        File::Spec->catfile($Bin, 'DiffGWASDeps', 'gnuplot', 'select_top_hits_from_wide.pl'),
        '--input', $args{wide_data},
        '--output', $out_tsv,
        '--focus-pvar', ($runner->{TOP_HIT_FOCUS_PVAR} || $runner->{MANHATTAN_P_VAR} || 'P'),
        '--focus-prefix', (($runner->{TOP_HIT_FOCUS_PVAR} || '') =~ /^(.*?)_(?:STD_DIFF_P|STD_P|DIFF_P|GROUP1_P|GROUP2_P)$/ ? $1 : ''),
        '--thresholds', $differential_thresholds,
        '--top-hit-dist-bp', ($runner->{TOP_HIT_DIST_BP} || '1e8'),
        '--max-hits', 15,
        '--maf-threshold', runner_maf_threshold($runner),
    );
    push @cmd, ('--gnomad-freq-file', $runner->{TOP_HIT_GNOMAD_FREQ_FILE})
        if defined $runner->{TOP_HIT_GNOMAD_FREQ_FILE} && length $runner->{TOP_HIT_GNOMAD_FREQ_FILE};
    push @cmd, ('--gnomad-pop-map', $runner->{TOP_HIT_GNOMAD_POP_MAP})
        if defined $runner->{TOP_HIT_GNOMAD_POP_MAP} && length $runner->{TOP_HIT_GNOMAD_POP_MAP};
    push @cmd, '--remove-x-chr' if $args{remove_x_chr};
    run_cmd(\@cmd, 'differential top-hit selection');
    return read_hits_tsv($out_tsv);
}

sub merge_hit_lists_for_gunplot {
    my (%args) = @_;
    my %seen;
    my @hits;
    for my $list (@{ $args{lists} || [] }) {
        for my $hit (@{ $list || [] }) {
            next unless ref($hit) eq 'HASH';
            my $snp = uc($hit->{SNP} || '');
            next unless length $snp;
            next if $seen{$snp}++;
            push @hits, { %{$hit} };
        }
    }
    if (defined($args{max_hits}) && $args{max_hits} =~ /^\d+$/ && $args{max_hits} > 0 && @hits > $args{max_hits}) {
        @hits = @hits[0 .. $args{max_hits} - 1];
    }
    for my $i (0 .. $#hits) {
        $hits[$i]{hit_order} = $i + 1;
    }
    return @hits;
}

sub resolve_ld_snps_for_query {
    my (%args) = @_;
    my %query = map { lc(trim($_)) => 1 } @{ $args{query_snps} || [] };
    my (%seen, %r2_for, @ld);
    my $add = sub {
        my ($snp, $r2) = @_;
        $snp = trim($snp // '');
        return unless length $snp;
        return if $query{lc $snp};
        my $key = lc $snp;
        if (defined($r2) && $r2 =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i) {
            $r2 = 0 + $r2;
            $r2_for{$key} = $r2 if !exists($r2_for{$key}) || $r2 > $r2_for{$key};
        }
        return if $seen{$key}++;
        push @ld, $snp;
    };
    for my $item (split /,/, ($args{explicit} // '')) {
        my ($snp, $r2) = split /:/, $item, 2;
        $add->($snp, $r2);
    }
    for my $item (split /,/, ($args{explicit_r2} // '')) {
        my ($snp, $r2) = split /:/, $item, 2;
        $add->($snp, $r2);
    }
    if (length(trim($args{explicit} // ''))) {
        print "[prep] User-supplied LD SNP marker overlay: " . join(',', @ld) . "\n";
        return (\@ld, \%r2_for);
    }

    my $runner = $args{runner} || {};
    my $population = uc($args{ld_population} || $runner->{LOCAL_LD_POPULATION} || 'EUR');
    my $min_r2 = 0 + ($args{ld_r2_threshold} || $runner->{LOCAL_LD_R2_THRESHOLD} || 0.8);
    my $audit_file = $args{audit_file} || '';
    $audit_file = localize_path($audit_file) if length $audit_file;
    if (!length($audit_file) && length($runner->{TOP_HIT_LD_AUDIT_BASENAME} || '')) {
        $audit_file = File::Spec->catfile($args{output_dir}, $runner->{TOP_HIT_LD_AUDIT_BASENAME});
    }
    if (length($audit_file) && -s $audit_file) {
        open my $fh, '<', $audit_file or die "Cannot read LD audit $audit_file: $!\n";
        my $header = <$fh> // '';
        chomp $header;
        $header =~ s/\r$//;
        my @cols = split /\t/, $header, -1;
        my %idx = map { lc(trim($cols[$_])) => $_ } 0 .. $#cols;
        if (exists $idx{candidate_snp} && exists $idx{lead_snp}) {
            while (my $line = <$fh>) {
                chomp $line;
                $line =~ s/\r$//;
                my @f = split /\t/, $line, -1;
                my $lead = $f[$idx{lead_snp}] // '';
                next unless $query{lc trim($lead)};
                if (exists $idx{selection_action}) {
                    my $action = $f[$idx{selection_action}] // '';
                    next unless $action =~ /^(?:PRUNED_LD|LD_LINKED)$/i;
                }
                # A clumping audit without r2 cannot establish that a proxy
                # meets the separate high-LD display threshold.
                next unless exists $idx{proxy_r2};
                my $r2 = $f[$idx{proxy_r2}] // '';
                next unless $r2 =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i
                    && $r2 >= $min_r2;
                $add->($f[$idx{candidate_snp}], $r2);
            }
        }
        close $fh;
    }

    my $cache_file = localize_path($args{ld_cache} || $runner->{LOCAL_LD_CACHE_TSV}
        || $runner->{TOP_HIT_LD_CACHE_TSV} || '');
    if (length($cache_file) && -s $cache_file) {
        open my $fh, '<', $cache_file or die "Cannot read LD cache $cache_file: $!\n";
        my $header = <$fh> // '';
        chomp $header;
        $header =~ s/\r$//;
        my @cols = split /\t/, $header, -1;
        my %idx = map { lc(trim($cols[$_])) => $_ } 0 .. $#cols;
        if (exists $idx{query_snp} && exists $idx{proxy_snp}) {
            while (my $line = <$fh>) {
                chomp $line;
                $line =~ s/\r$//;
                my @f = split /\t/, $line, -1;
                my $lead = $f[$idx{query_snp}] // '';
                next unless $query{lc trim($lead)};
                if (exists $idx{ld_population}) {
                    next unless uc(trim($f[$idx{ld_population}] // '')) eq $population;
                }
                if (exists $idx{proxy_r2}) {
                    my $r2 = $f[$idx{proxy_r2}] // '';
                    next unless $r2 =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i && $r2 >= $min_r2;
                }
                $add->($f[$idx{proxy_snp}], (exists($idx{proxy_r2}) ? $f[$idx{proxy_r2}] : undef));
            }
        }
        close $fh;
    }

    my $helper = File::Spec->catfile($Bin, 'DiffGWASDeps', 'resolve_haploreg_high_ld.pl');
    if (-f $helper) {
        my $web_cache = File::Spec->catfile(
            $Bin, 'cache', 'haploreg_ld', 'high_ld_queries.tsv'
        );
        my @cmd = (
            $^X, $helper,
            '--query-snps', join(',', @{ $args{query_snps} || [] }),
            '--population', $population,
            '--min-r2', $min_r2,
            '--web-cache', $web_cache,
        );
        push @cmd, ('--local-cache', $cache_file) if length($cache_file) && -s $cache_file;
        push @cmd, '--no-web-fallback' unless $args{ld_web_fallback};
        if (open my $pipe, '-|', @cmd) {
            while (my $line = <$pipe>) {
                chomp $line;
                if ($line =~ /^LD_SNPS\t(.*)$/) {
                    $add->($_) for split /,/, $1;
                }
                elsif ($line =~ /^LD_R2_PAIRS\t(.*)$/) {
                    for my $item (split /,/, $1) {
                        my ($snp, $r2) = split /:/, $item, 2;
                        $add->($snp, $r2);
                    }
                }
            }
            my $ok = close $pipe;
            warn "WARNING: High-LD resolver exited unsuccessfully; continuing without automatic stars.\n"
                unless $ok;
        }
        else {
            warn "WARNING: Cannot start high-LD resolver $helper: $!\n";
        }
    }
    print "[prep] LD-linked SNP marker overlay for " . join(',', @{ $args{query_snps} || [] })
        . " ($population, r2 >= $min_r2): " . join(',', @ld) . "\n" if @ld;
    print "[warn] No high-LD proxies resolved for " . join(',', @{ $args{query_snps} || [] })
        . " ($population, r2 >= $min_r2); no LD marker overlay will be drawn.\n" unless @ld;
    return (\@ld, \%r2_for);
}

sub plot_local_series {
    my (%args) = @_;
    my $runner = $args{runner};
    my $base_name = $args{output_base};
    $base_name =~ s/\.html$//i;

    my @images;
    my @parts;
    my $batch_size = $args{batch_size} || scalar(@{ $args{hits} }) || 1;
    my $panel_columns = $args{panel_columns} || 0;
    my $cache_dir = $args{gtf_cache_dir}
        || File::Spec->catfile($args{output_dir}, '.gnuplot_gtf_cache');
    if (!$args{gtf_cache_dir} && !-e $cache_dir) {
        my $legacy_cache_dir = File::Spec->catfile($args{output_dir}, '.gunplot_gtf_cache');
        if (-d $legacy_cache_dir) {
            if (rename($legacy_cache_dir, $cache_dir)) {
                print "[migrate] Renamed legacy GTF cache directory $legacy_cache_dir to $cache_dir\n";
            }
            else {
                warn "WARNING: Could not rename legacy GTF cache directory $legacy_cache_dir to $cache_dir: $!. A new canonical cache will be created.\n";
            }
        }
    }
    make_path($cache_dir) unless -d $cache_dir;
    # Resolve each target from the indexed GWAS source before grouping.  On a
    # cold run the caller normally supplies only SNP names, so grouping before
    # this step silently treated nearby targets as unrelated loci.
    my $locus_sources = prepare_locus_wide_sources(
        hits          => $args{hits},
        output_dir    => $args{output_dir},
        window_bp     => $args{window_bp},
        source_long   => $args{source_long},
        preset_config => $args{preset_config},
        force         => $args{force},
    );

    my (%label_snps_for_index, %skip_hit_index);
    if ($args{kind} eq 'local_gtf' || $args{kind} eq 'local_manhattan') {
        my $window = 0 + $args{window_bp};
        for my $i (0 .. $#{ $args{hits} }) {
            next if $skip_hit_index{$i};
            my $anchor = $args{hits}[$i];
            my @members = ($anchor->{SNP});
            if (defined($anchor->{CHR}) && defined($anchor->{BP})) {
                for my $j ($i + 1 .. $#{ $args{hits} }) {
                    next if $skip_hit_index{$j};
                    my $candidate = $args{hits}[$j];
                    next unless defined($candidate->{CHR}) && defined($candidate->{BP});
                    next unless normalize_chr_label($candidate->{CHR}) eq normalize_chr_label($anchor->{CHR});
                    next unless abs((0 + $candidate->{BP}) - (0 + $anchor->{BP})) <= $window;
                    push @members, $candidate->{SNP};
                    $skip_hit_index{$j} = 1;
                }
            }
            $label_snps_for_index{$i} = \@members;
            if (@members > 1) {
                print "[prep] Combining close $args{kind} query SNPs in one shared window centered on $anchor->{SNP}: "
                    . join(', ', @members) . "\n";
            }
        }
    }
    my $n_hits = scalar grep { !$skip_hit_index{$_} }
        0 .. $#{ $args{hits} };
    if ($args{kind} eq 'local_gtf' || $args{kind} eq 'local_manhattan') {
        for my $i (keys %label_snps_for_index) {
            my $members = $label_snps_for_index{$i};
            next unless @{$members} > 1;
            my $shared = $locus_sources->{uc($members->[0])};
            next unless $shared;
            $locus_sources->{uc($_)} = $shared for @{$members}[1 .. $#{$members}];
        }
    }

    write_top_hits_csv(
        path       => File::Spec->catfile($args{output_dir}, $args{top_csv_name}),
        hits       => $args{hits},
        batch_size => $batch_size,
        runner     => $runner,
        wide_data  => $args{wide_data},
        locus_sources => $locus_sources,
    );

    my $render_idx = -1;
    for my $hit_idx (0 .. $#{ $args{hits} }) {
        next if $skip_hit_index{$hit_idx};
        $render_idx++;
        my $hit = $args{hits}[$hit_idx];
        my @label_snps = @{ $label_snps_for_index{$hit_idx} || [$hit->{SNP}] };
        my $label_snps_csv = join(',', @label_snps);
        my $locus_title_snps = join(', ', @label_snps);
        my ($ld_snps_ref, $ld_r2_ref) = $args{highlight_high_ld_snps} ? resolve_ld_snps_for_query(
            explicit   => $args{ld_snps},
            explicit_r2=> $args{ld_r2_values},
            audit_file => $args{ld_audit_file},
            ld_cache   => $args{ld_cache},
            ld_population => $args{ld_population},
            ld_r2_threshold => $args{ld_r2_threshold},
            ld_web_fallback => $args{ld_web_fallback},
            runner     => $runner,
            output_dir => $args{output_dir},
            query_snps => \@label_snps,
        ) : ([], {});
        my @ld_snps = @{ $ld_snps_ref || [] };
        my $ld_snps_csv = join(',', @ld_snps);
        my $ld_r2_values = join(',', map {
            my $key = lc $_;
            exists($ld_r2_ref->{$key}) ? $_ . ':' . sprintf('%.6g', $ld_r2_ref->{$key}) : ()
        } @ld_snps);
        my $safe_snp = safe_name($hit->{SNP});
        my $safe_window = safe_name($args{window_bp});
        my $batch_index = int($render_idx / $batch_size);
        my $batch_start = $batch_index * $batch_size;
        my $batch_end = $batch_start + $batch_size - 1;
        $batch_end = $n_hits - 1 if $batch_end > $n_hits - 1;
        my $batch_count = $batch_end - $batch_start + 1;
        my $batch_cols = compute_panel_columns($batch_count, $panel_columns);
        my $batch_annotation_mode = normalize_local_manhattan_annotation_mode($args{annotation_mode}, $batch_count);
        my $batch_pos = $render_idx - $batch_start;
        my $batch_col = $batch_pos % $batch_cols;
        my $locus_prefix = File::Spec->catfile($args{output_dir}, $base_name . '_' . $safe_snp);
        my $existing_locus_manifest = $locus_prefix . '.manifest.tsv';
        if (-s $existing_locus_manifest && (!defined $hit->{CHR} || !defined $hit->{BP})) {
            my $prior_metrics = read_manifest_tsv($existing_locus_manifest);
            my $prior_snp = $prior_metrics->{snp} // '';
            if (lc($prior_snp) eq lc($hit->{SNP})) {
                $hit->{CHR} = $prior_metrics->{chr} if !defined($hit->{CHR}) && defined($prior_metrics->{chr});
                $hit->{BP}  = $prior_metrics->{bp}  if !defined($hit->{BP})  && defined($prior_metrics->{bp});
            }
        }
        my @required_locus_outputs = (
            $locus_prefix . '.png',
            $locus_prefix . '.plot.tsv',
            $locus_prefix . '.manifest.tsv',
        );
        if ($args{with_gtf} || $batch_annotation_mode eq 'gtf') {
            push @required_locus_outputs, $locus_prefix . '.genes.tsv';
        }
        my $expected_gtf_file = infer_cached_locus_gtf_path(
            output_dir => $args{output_dir},
            snp        => $hit->{SNP},
            window_bp  => $args{window_bp},
            runner     => $runner,
        );
        my $expected_has_gtf = ($args{with_gtf} || $batch_annotation_mode eq 'gtf') ? 1 : 0;
        my $expected_bottom_gene = '';
        if ($args{kind} eq 'local_manhattan' && !$expected_has_gtf) {
            $expected_bottom_gene = $hit->{gene} || '';
            $expected_bottom_gene = extract_gene_from_snp_gene($hit->{snp_gene}) if !length($expected_bottom_gene);
            $expected_bottom_gene = '' if $expected_bottom_gene =~ /^(?:NA|N\/A|null|\.)$/i;
        }
        my ($cache_reusable, $cache_reason) = local_locus_cache_is_reusable(
            required          => \@required_locus_outputs,
            manifest          => $existing_locus_manifest,
            png               => $locus_prefix . '.png',
            renderer          => File::Spec->catfile($Bin, 'DiffGWASDeps', 'gnuplot', 'pdl_gunplot_local_locus.pl'),
            snp               => $hit->{SNP},
            window_bp         => $args{window_bp},
            pcols             => join(',', @{ $args{pcols} }),
            zcols             => join(',', @{ $args{zcols} || [] }),
            labels            => join('|', @{ $args{labels} }),
            title             => sprintf('%s: %s (%s:%s)', $args{html_title}, $locus_title_snps, ($hit->{CHR} // ''), ($hit->{BP} // '')),
            label_snps        => $label_snps_csv,
            ld_snps           => $ld_snps_csv,
            ld_marker_symbol  => $args{ld_marker_symbol},
            ld_marker_color   => $args{ld_marker_color},
            ld_display_mode   => $args{ld_display_mode},
            ld_r2_values      => $ld_r2_values,
            ld_heatmap_colors => $args{ld_heatmap_colors},
            ld_population     => $args{ld_population},
            height            => $args{height},
            sig               => ($runner->{TOP_HIT_SIGNAL_THRSHD} || '1e-6'),
            has_gtf           => $expected_has_gtf,
            gtf_file          => $expected_gtf_file,
            hide_y_axis       => ($args{kind} eq 'local_manhattan' && $batch_col > 0) ? 1 : 0,
            bottom_snp_label  => ($args{kind} eq 'local_manhattan' && !$expected_has_gtf) ? $hit->{SNP} : '',
            bottom_gene_label => $expected_bottom_gene,
        );
        my $need_locus_render = $args{force} || !$cache_reusable;
        if ($need_locus_render && !$args{force} && all_nonempty_files(@required_locus_outputs)) {
            print "[cache] regenerating gunplot $args{kind} locus artifacts for $hit->{SNP}: $cache_reason\n";
        }

        if (!$need_locus_render) {
            my $manifest_metrics = read_manifest_tsv($locus_prefix . '.manifest.tsv');
            $hit->{CHR} = $manifest_metrics->{TARGET_CHR} || $manifest_metrics->{target_chr} || $manifest_metrics->{CHR} || $manifest_metrics->{chr}
                if ($manifest_metrics->{TARGET_CHR} || $manifest_metrics->{target_chr} || $manifest_metrics->{CHR} || $manifest_metrics->{chr});
            $hit->{BP}  = $manifest_metrics->{TARGET_BP} || $manifest_metrics->{target_bp} || $manifest_metrics->{BP} || $manifest_metrics->{bp}
                if ($manifest_metrics->{TARGET_BP} || $manifest_metrics->{target_bp} || $manifest_metrics->{BP} || $manifest_metrics->{bp});
            ensure_genes_tsv_from_gtf(
                gtf_file  => $expected_gtf_file,
                gene_tsv  => $locus_prefix . '.genes.tsv',
            ) if $expected_gtf_file;
            print "[skip] reusing existing gunplot $args{kind} locus artifacts for $hit->{SNP}\n";
            push @images, {
                slot_index=> $batch_pos,
                snp      => $locus_title_snps,
                chr      => $hit->{CHR},
                bp       => $hit->{BP},
                gene     => ($hit->{gene} || ''),
                snp_gene => ($hit->{snp_gene} || ''),
                image    => $locus_prefix . '.png',
                plot_tsv => $locus_prefix . '.plot.tsv',
                gtf_file => $expected_gtf_file,
                manifest => $locus_prefix . '.manifest.tsv',
            };
            next;
        }

        my $locus_input = $args{wide_data};
        my $cached_locus_input = File::Spec->catfile($args{output_dir}, "gunplot_locus_${safe_snp}_window_${safe_window}.wide.tsv.gz");
        my $cached_locus_manifest = File::Spec->catfile($args{output_dir}, "gunplot_locus_${safe_snp}_window_${safe_window}.wide.manifest.tsv");
        my $legacy_locus_input = File::Spec->catfile($args{output_dir}, $base_name . '_' . $safe_snp . '.wide.tsv.gz');
        my $legacy_locus_manifest = File::Spec->catfile($args{output_dir}, $base_name . '_' . $safe_snp . '.wide.manifest.tsv');
        if (!(-s $cached_locus_input && -s $cached_locus_manifest) && -s $legacy_locus_input && -s $legacy_locus_manifest) {
            $cached_locus_input = $legacy_locus_input;
            $cached_locus_manifest = $legacy_locus_manifest;
        }
        if (!(-s $cached_locus_input && -s $cached_locus_manifest)) {
            for my $candidate (glob(File::Spec->catfile($args{output_dir}, '*' . $safe_snp . '.wide.tsv.gz'))) {
                (my $candidate_manifest = $candidate) =~ s/\.wide\.tsv\.gz$/.wide.manifest.tsv/;
                next unless -s $candidate && -s $candidate_manifest;
                $cached_locus_input = $candidate;
                $cached_locus_manifest = $candidate_manifest;
                last;
            }
        }
        if (-s $cached_locus_input && -s $cached_locus_manifest) {
            $locus_input = $cached_locus_input;
            my $manifest_metrics = read_manifest_tsv($cached_locus_manifest);
            $hit->{CHR} = $manifest_metrics->{TARGET_CHR} || $manifest_metrics->{target_chr} || $manifest_metrics->{CHR} || $manifest_metrics->{chr}
                if ($manifest_metrics->{TARGET_CHR} || $manifest_metrics->{target_chr} || $manifest_metrics->{CHR} || $manifest_metrics->{chr});
            $hit->{BP}  = $manifest_metrics->{TARGET_BP} || $manifest_metrics->{target_bp} || $manifest_metrics->{BP} || $manifest_metrics->{bp}
                if ($manifest_metrics->{TARGET_BP} || $manifest_metrics->{target_bp} || $manifest_metrics->{BP} || $manifest_metrics->{bp});
        }
        if ((!defined $hit->{CHR} || !defined $hit->{BP}) && $args{wide_data}) {
            my ($chr, $bp) = locate_snp_in_gzip($args{wide_data}, $hit->{SNP});
            if (defined $chr && defined $bp) {
                $hit->{CHR} = $chr;
                $hit->{BP} = $bp;
            }
        }
        if ($args{source_long} && -s $args{source_long} && $args{preset_config} && -f $args{preset_config} && !(-s $cached_locus_input && -s $cached_locus_manifest)) {
            $locus_input = $cached_locus_input;
            my $locus_manifest = $cached_locus_manifest;
            my @extract_cmd = (
                $^X,
                File::Spec->catfile($Bin, 'DiffGWASDeps', 'extract_single_snp_wide_diff_gwas.pl'),
                '--config', $args{preset_config},
                '--input', $args{source_long},
                '--target-snp', $hit->{SNP},
                '--window-bp', $args{window_bp},
                '--output', $locus_input,
                '--manifest', $locus_manifest,
                '--output-dir', $args{output_dir},
            );
            push @extract_cmd, ('--target-chr', $hit->{CHR}) if defined $hit->{CHR} && $hit->{CHR} ne '';
            push @extract_cmd, ('--target-bp',  $hit->{BP})  if defined $hit->{BP}  && $hit->{BP}  ne '';
            run_cmd(\@extract_cmd, "local wide extraction for $hit->{SNP}");
            my $manifest_metrics = read_manifest_tsv($locus_manifest);
            $hit->{CHR} = $manifest_metrics->{TARGET_CHR} || $manifest_metrics->{target_chr} || $manifest_metrics->{CHR} || $manifest_metrics->{chr}
                if ($manifest_metrics->{TARGET_CHR} || $manifest_metrics->{target_chr} || $manifest_metrics->{CHR} || $manifest_metrics->{chr});
            $hit->{BP}  = $manifest_metrics->{TARGET_BP} || $manifest_metrics->{target_bp} || $manifest_metrics->{BP} || $manifest_metrics->{bp}
                if ($manifest_metrics->{TARGET_BP} || $manifest_metrics->{target_bp} || $manifest_metrics->{BP} || $manifest_metrics->{bp});
        }
        if (-s $cached_locus_input && -s $cached_locus_manifest) {
            $locus_input = $cached_locus_input;
        }
        my $gtf_file = '';
        if ($args{with_gtf} || $args{kind} eq 'local_manhattan') {
            die "CHR/BP metadata are missing for $hit->{SNP}; locus extraction did not resolve coordinates.\n"
                unless defined $hit->{CHR} && defined $hit->{BP};
            my $npc_flag = $runner->{LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES} ? 1 : 0;
            $gtf_file = File::Spec->catfile(
                $args{output_dir},
                "gunplot_locus_${safe_snp}_window_${safe_window}_npc${npc_flag}.gtf.tsv"
            );
            my $region_start = $hit->{BP} - (0 + $args{window_bp});
            $region_start = 1 if $region_start < 1;
            my $region_end = $hit->{BP} + (0 + $args{window_bp});
            if (!$args{force} && -s $gtf_file) {
                print "[skip] reusing cached gunplot GTF subset $gtf_file\n";
            }
            else {
                my @gtf_cmd = (
                    $^X,
                    File::Spec->catfile($Bin, 'DiffGWASDeps', 'extract_gencode_gtf_subset.pl'),
                    '--cache-dir', $cache_dir,
                    '--output', $gtf_file,
                    '--region', "$hit->{CHR}:$region_start:$region_end",
                );
                push @gtf_cmd, ('--reference-build', $runner->{REFERENCE_BUILD})
                    if defined $runner->{REFERENCE_BUILD} && length $runner->{REFERENCE_BUILD};
                push @gtf_cmd, ('--gtf-url', $runner->{GTF_GZ_URL})
                    if defined $runner->{GTF_GZ_URL} && length $runner->{GTF_GZ_URL};
                if (!$runner->{LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}) {
                    push @gtf_cmd, '--no-include-non-protein-coding';
                }
                run_cmd(\@gtf_cmd, "GTF subset for $hit->{SNP}");
            }
        }

        my @cmd = (
            $^X,
            File::Spec->catfile($Bin, 'DiffGWASDeps', 'gnuplot', 'pdl_gunplot_local_locus.pl'),
            '--data', $locus_input,
            '--snp', $hit->{SNP},
            '--label-snps', $label_snps_csv,
            '--out-prefix', $locus_prefix,
            '--window-bp', $args{window_bp},
            '--pcols', join(',', @{ $args{pcols} }),
            '--labels', join('|', @{ $args{labels} }),
            '--title', sprintf('%s: %s (%s:%s)', $args{html_title}, $locus_title_snps, $hit->{CHR}, $hit->{BP}),
            '--height', $args{height},
            '--gnuplot', $args{gnuplot},
            '--sig', ($runner->{TOP_HIT_SIGNAL_THRSHD} || '1e-6'),
        );
        push @cmd, ('--ld-snps', $ld_snps_csv) if length $ld_snps_csv;
        if (length $ld_snps_csv) {
            push @cmd, ('--ld-marker-symbol', ($args{ld_marker_symbol} || 'star'));
            push @cmd, ('--ld-marker-color', ($args{ld_marker_color} || 'black'));
            push @cmd, ('--ld-display-mode', ($args{ld_display_mode} || 'markers'));
            push @cmd, ('--ld-population', ($args{ld_population} || 'EUR'));
            push @cmd, ('--ld-heatmap-colors', ($args{ld_heatmap_colors} || '#f7fbff,#6baed6,#54278f'));
            push @cmd, ('--ld-r2-values', $ld_r2_values) if length $ld_r2_values;
        }
        if ($args{kind} eq 'local_manhattan') {
            my $bottom_gene = $hit->{gene} || '';
            $bottom_gene = extract_gene_from_snp_gene($hit->{snp_gene}) if !length($bottom_gene);
            $bottom_gene = '' if $bottom_gene =~ /^(?:NA|N\/A|null|\.)$/i;
            if ($batch_annotation_mode eq 'gtf') {
                push @cmd, ('--gtf', $gtf_file) if $gtf_file;
                push @cmd, '--hide-y-axis' if $batch_col > 0;
            }
            else {
                push @cmd, ('--bottom-snp-label', $hit->{SNP});
                push @cmd, ('--bottom-gene-label', $bottom_gene) if length $bottom_gene;
                push @cmd, '--hide-y-axis' if $batch_col > 0;
            }
        }
        if ($args{zcols} && ref($args{zcols}) eq 'ARRAY' && @{ $args{zcols} }) {
            push @cmd, ('--zcols', join(',', @{ $args{zcols} }));
        }
        push @cmd, ('--gtf', $gtf_file) if $args{with_gtf} && $gtf_file;
        run_cmd(\@cmd, "$args{kind} locus plot for $hit->{SNP}");
        ensure_genes_tsv_from_gtf(
            gtf_file => $gtf_file,
            gene_tsv => $locus_prefix . '.genes.tsv',
        ) if $gtf_file;

        push @images, {
            slot_index=> $batch_pos,
            snp      => $locus_title_snps,
            chr      => $hit->{CHR},
            bp       => $hit->{BP},
            gene     => ($hit->{gene} || ''),
            snp_gene => ($hit->{snp_gene} || ''),
            image    => $locus_prefix . '.png',
            plot_tsv => $locus_prefix . '.plot.tsv',
            gtf_file => $gtf_file,
            manifest => $locus_prefix . '.manifest.tsv',
        };
    }

    my @batches;
    while (@images) {
        push @batches, [ splice(@images, 0, $batch_size) ];
    }

    my $main_html = File::Spec->catfile($args{output_dir}, $base_name . '.html');
    if ($args{kind} eq 'local_manhattan') {
        my @batch_images;
        for my $i (0 .. $#batches) {
            my $batch_png = File::Spec->catfile(
                $args{output_dir},
                $i == 0 ? ($base_name . '.png') : sprintf('%s_part%d.png', $base_name, $i + 1)
            );
            my $batch_cols = compute_panel_columns(scalar(@{ $batches[$i] }), $panel_columns);
            my $batch_annotation_mode = normalize_local_manhattan_annotation_mode($args{annotation_mode}, scalar(@{ $batches[$i] }));
            # A group of adjacent query SNPs is already rendered as one
            # coordinate-aware locus by pdl_gunplot_local_locus.pl.  Rebuilding
            # that single locus through the multi-locus compositor discarded
            # its query-SNP labels and made old two-panel output easy to reuse.
            # Publish the authoritative locus image directly under the stable
            # combined-output name instead.
            if (@{ $batches[$i] } == 1) {
                my $item = $batches[$i][0];
                my $manifest = read_manifest_tsv($item->{manifest});
                my $label_snps = $manifest->{label_snps} || '';
                if ($label_snps =~ /,/) {
                    if ($args{force}
                        || !-s $batch_png
                        || !target_is_newer_than_inputs($batch_png, $item->{image}, $item->{manifest})) {
                        copy($item->{image}, $batch_png)
                            or die "Cannot publish grouped local Manhattan image $item->{image} as $batch_png: $!\n";
                        print "[done] published one grouped local Manhattan locus with labels $label_snps: $batch_png\n";
                    }
                    else {
                        print "[skip] reusing grouped single-locus local Manhattan batch $batch_png\n";
                    }
                    push @batch_images, $batch_png;
                    next;
                }
            }
            if ($batch_annotation_mode eq 'gtf') {
                if (!$args{force} && -s $batch_png && target_is_newer_than_inputs($batch_png, batch_dependency_files($batches[$i], 1))) {
                    print "[skip] reusing existing combined gunplot local Manhattan GTF batch $batch_png\n";
                }
                else {
                    render_combined_local_manhattan_gtf_batch(
                        output_png  => $batch_png,
                        output_dir  => $args{output_dir},
                        output_base => ($i == 0 ? $base_name : sprintf('%s_part%d', $base_name, $i + 1)),
                        title       => sprintf('%s%s', $args{html_title}, @batches > 1 ? sprintf(' (Part %d)', $i + 1) : ''),
                        images      => $batches[$i],
                        labels      => $args{labels},
                        columns     => $batch_cols,
                        top_logp    => 8,
                        sig_y       => safe_neglog10_text($runner->{TOP_HIT_SIGNAL_THRSHD} || '1e-6'),
                        gnuplot     => $args{gnuplot},
                    );
                }
            } else {
                if (!$args{force} && -s $batch_png && target_is_newer_than_inputs($batch_png, batch_dependency_files($batches[$i], 0))) {
                    print "[skip] reusing existing combined gunplot local Manhattan batch $batch_png\n";
                }
                else {
                    render_combined_local_manhattan_batch(
                        output_png  => $batch_png,
                        output_dir  => $args{output_dir},
                        output_base => ($i == 0 ? $base_name : sprintf('%s_part%d', $base_name, $i + 1)),
                        title       => sprintf('%s%s', $args{html_title}, @batches > 1 ? sprintf(' (Part %d)', $i + 1) : ''),
                        images      => $batches[$i],
                        labels      => $args{labels},
                        columns     => $batch_cols,
                        annotation_mode => $args{annotation_mode},
                        top_logp    => 8,
                        sig_y       => safe_neglog10_text($runner->{TOP_HIT_SIGNAL_THRSHD} || '1e-6'),
                        gnuplot     => $args{gnuplot},
                    );
                }
            }
            push @batch_images, $batch_png;
        }

        if (@batch_images == 1) {
            write_single_image_html(
                path      => $main_html,
                title     => $args{html_title},
                image_rel => basename($batch_images[0]),
            );
        }
        else {
            for my $i (0 .. $#batch_images) {
                my $part_html = File::Spec->catfile($args{output_dir}, sprintf('%s_part%d.html', $base_name, $i + 1));
                write_single_image_html(
                    path      => $part_html,
                    title     => sprintf('%s (Part %d)', $args{html_title}, $i + 1),
                    image_rel => basename($batch_images[$i]),
                );
                push @parts, $part_html;
            }
            write_index_html(
                path      => $main_html,
                title     => $args{html_title},
                parts     => \@parts,
                top_csv   => gunplotize_name($runner->{LOCAL_TOP_HITS_CSV_BASENAME} || 'gunplot_top_hits.csv'),
                batch_msg => sprintf('Top hits were split into %d combined gunplot batch figures.', scalar(@parts)),
            );
        }
    }
    else {
        if (@batches == 1) {
            write_gallery_html(
                path    => $main_html,
                title   => $args{html_title},
                images  => $batches[0],
                top_csv => gunplotize_name($runner->{LOCAL_TOP_HITS_CSV_BASENAME} || 'gunplot_top_hits.csv'),
            );
        }
        else {
            for my $i (0 .. $#batches) {
                my $part_html = File::Spec->catfile($args{output_dir}, sprintf('%s_part%d.html', $base_name, $i + 1));
                write_gallery_html(
                    path    => $part_html,
                    title   => sprintf('%s (Part %d)', $args{html_title}, $i + 1),
                    images  => $batches[$i],
                    top_csv => gunplotize_name($runner->{LOCAL_TOP_HITS_CSV_BASENAME} || 'gunplot_top_hits.csv'),
                );
                push @parts, $part_html;
            }
            write_index_html(
                path      => $main_html,
                title     => $args{html_title},
                parts     => \@parts,
                top_csv   => gunplotize_name($runner->{LOCAL_TOP_HITS_CSV_BASENAME} || 'gunplot_top_hits.csv'),
                batch_msg => sprintf('Top hits were split into %d gunplot batch pages.', scalar(@parts)),
            );
        }
    }

    return (
        $args{kind} . '_html' => $main_html,
        ($args{kind} eq 'local_manhattan'
            ? ($args{kind} . '_png' => File::Spec->catfile($args{output_dir}, $base_name . '.png'))
            : ()),
    );
}

sub compute_panel_columns {
    my ($n_images, $requested) = @_;
    $n_images ||= 1;
    if ($requested && $requested =~ /^\d+$/ && $requested > 0) {
        return $requested > $n_images ? $n_images : $requested;
    }
    my $cols = int(sqrt($n_images));
    $cols++ if $cols * $cols < $n_images;
    $cols = 1 if $cols < 1;
    $cols = $n_images if $cols > $n_images;
    return $cols;
}

sub all_nonempty_files {
    for my $path (@_) {
        return 0 unless defined $path && -s $path;
    }
    return 1;
}

sub local_locus_cache_is_reusable {
    my (%args) = @_;
    return (0, 'one or more expected output files are absent or empty')
        unless all_nonempty_files(@{ $args{required} || [] });

    my $manifest = $args{manifest};
    return (0, 'plot manifest is absent or empty') unless defined $manifest && -s $manifest;
    my $metrics = read_manifest_tsv($manifest);
    return (0, 'legacy plot manifest has no cache schema')
        unless defined $metrics->{cache_schema} && $metrics->{cache_schema} =~ /^\d+$/ && $metrics->{cache_schema} >= 2;

    my @checks = (
        ['snp',               ($args{snp} // '')],
        ['pcols',             ($args{pcols} // '')],
        ['label_snps',        ($args{label_snps} // $args{snp} // '')],
        ['ld_snps',           ($args{ld_snps} // '')],
        ['ld_marker_symbol',  ($args{ld_marker_symbol} // 'star')],
        ['ld_marker_color',   ($args{ld_marker_color} // 'black')],
        ['ld_display_mode',   ($args{ld_display_mode} // 'none')],
        ['ld_r2_values',      ($args{ld_r2_values} // '')],
        ['ld_heatmap_colors', ($args{ld_heatmap_colors} // '#f7fbff,#6baed6,#54278f')],
        ['ld_population',     ($args{ld_population} // 'EUR')],
        ['zcols',             ($args{zcols} // '')],
        ['labels',            ($args{labels} // '')],
        ['title',             ($args{title} // '')],
        ['height',            ($args{height} // '')],
        ['sig',               ($args{sig} // '')],
        ['has_gtf',           ($args{has_gtf} ? 1 : 0)],
        ['hide_y_axis',       ($args{hide_y_axis} ? 1 : 0)],
        ['bottom_snp_label',  ($args{bottom_snp_label} // '')],
        ['bottom_gene_label', ($args{bottom_gene_label} // '')],
    );
    for my $check (@checks) {
        my ($key, $expected) = @$check;
        my $actual = defined $metrics->{$key} ? $metrics->{$key} : '';
        return (0, "plot setting '$key' changed") unless $actual eq "$expected";
    }

    my $actual_window = $metrics->{window_bp};
    return (0, 'plot window is missing from manifest') unless defined $actual_window && $actual_window ne '';
    return (0, 'plot window changed') unless (0 + $actual_window) == (0 + $args{window_bp});

    if ($args{has_gtf}) {
        return (0, 'current target/window GTF subset is absent') unless $args{gtf_file} && -s $args{gtf_file};
        return (0, 'GTF subset path changed') unless cache_paths_equal($metrics->{gtf}, $args{gtf_file});
    }
    elsif (defined $metrics->{gtf} && length $metrics->{gtf}) {
        return (0, 'cached plot unexpectedly contains GTF annotation');
    }

    my @dependencies = grep { defined && length } (
        $metrics->{input},
        ($args{has_gtf} ? $args{gtf_file} : ''),
        $args{renderer},
    );
    return (0, 'a plot input is absent or newer than the cached PNG')
        unless target_is_newer_than_inputs($args{png}, @dependencies);
    return (1, 'cache provenance and dependencies match');
}

sub cache_paths_equal {
    my ($left, $right) = @_;
    return 0 unless defined $left && defined $right && length $left && length $right;
    for ($left, $right) {
        s{\\}{/}g;
        s{/+}{/}g;
        s{/$}{};
        $_ = lc($_) if $^O =~ /(?:cygwin|MSWin32)/i || /^[A-Za-z]:/;
    }
    return $left eq $right;
}

sub target_is_newer_than_inputs {
    my ($target, @inputs) = @_;
    return 0 unless defined $target && -s $target;
    my $target_mtime = (stat($target))[9];
    return 0 unless defined $target_mtime;
    for my $path (@inputs) {
        next unless defined $path && length $path;
        return 0 unless -e $path;
        my $mtime = (stat($path))[9];
        return 0 unless defined $mtime;
        return 0 if $mtime > $target_mtime;
    }
    return 1;
}

sub batch_dependency_files {
    my ($items, $with_gtf) = @_;
    my @deps;
    for my $item (@{ $items || [] }) {
        push @deps, grep { defined && length } (
            $item->{plot_tsv},
            $item->{manifest},
        );
    }
    return @deps;
}

sub compose_png_grid {
    my (%args) = @_;
    my @images = @{ $args{images} || [] };
    die "No images were supplied for compose_png_grid\n" unless @images;
    my $cols = $args{columns} || compute_panel_columns(scalar(@images), 0);
    $cols = 1 if $cols < 1;
    my $rows = int((scalar(@images) + $cols - 1) / $cols);
    my $pad = 18;

    my @gd_images;
    my ($max_w, $max_h) = (0, 0);
    for my $path (@images) {
        open my $fh, '<', $path or die "Cannot read $path: $!\n";
        binmode $fh;
        my $img = GD::Image->newFromPng($fh, 1)
            or die "Cannot decode PNG $path\n";
        close $fh;
        push @gd_images, $img;
        $max_w = $img->width  if $img->width  > $max_w;
        $max_h = $img->height if $img->height > $max_h;
    }

    my $canvas_w = $cols * $max_w + ($cols + 1) * $pad;
    my $canvas_h = $rows * $max_h + ($rows + 1) * $pad;
    my $canvas = GD::Image->newTrueColor($canvas_w, $canvas_h);
    my $white = $canvas->colorAllocate(255, 255, 255);
    $canvas->filledRectangle(0, 0, $canvas_w - 1, $canvas_h - 1, $white);
    $canvas->alphaBlending(1);
    $canvas->saveAlpha(1);

    for my $i (0 .. $#gd_images) {
        my $row = int($i / $cols);
        my $col = $i % $cols;
        my $x = $pad + $col * ($max_w + $pad);
        $x += int(($max_w - $gd_images[$i]->width) / 2);
        my $y = $pad + $row * ($max_h + $pad);
        $y += int(($max_h - $gd_images[$i]->height) / 2);
        $canvas->copy($gd_images[$i], $x, $y, 0, 0, $gd_images[$i]->width, $gd_images[$i]->height);
    }

    open my $out, '>', $args{output_png} or die "Cannot write $args{output_png}: $!\n";
    binmode $out;
    print {$out} $canvas->png;
    close $out or die "Cannot close $args{output_png}: $!\n";
}

sub render_combined_local_manhattan_batch {
    my (%args) = @_;
    my @items = @{ $args{images} || [] };
    die "No local Manhattan loci supplied for batch render\n" unless @items;
    my $n = scalar @items;
    my $cols = $args{columns} || $n;
    $cols = $n if $cols > $n;
    $cols = 1 if $cols < 1;
    my $annotation_mode = normalize_local_manhattan_annotation_mode($args{annotation_mode}, $n);

    my $combined_tsv = File::Spec->catfile($args{output_dir}, $args{output_base} . '.combined.tsv');
    my $gp_file = File::Spec->catfile($args{output_dir}, $args{output_base} . '.combined.gp');

    open my $out, '>', $combined_tsv or die "Cannot write $combined_tsv: $!\n";
    print {$out} join("\t", qw(X Y LOCUS)), "\n";

    for my $i (0 .. $#items) {
        my $plot_tsv = $items[$i]{plot_tsv};
        my $manifest = read_manifest_tsv($items[$i]{manifest});
        my $target_bp = $manifest->{BP} || $manifest->{bp} || $manifest->{TARGET_BP} || $manifest->{target_bp} || $items[$i]{bp};
        my $window_bp = $manifest->{WINDOW_BP} || $manifest->{window_bp} || 0;
        my $start = $manifest->{WINDOW_START} || $manifest->{window_start} || ($target_bp - $window_bp);
        my $end   = $manifest->{WINDOW_END} || $manifest->{window_end} || ($target_bp + $window_bp);
        $start = 1 if !defined $start || $start < 1;
        $end = $start + 1 if !defined $end || $end <= $start;
        open my $fh, '<', $plot_tsv or die "Cannot read $plot_tsv: $!\n";
        my $header = <$fh>;
        while (my $line = <$fh>) {
            chomp $line;
            $line =~ s/\r$//;
            next unless length $line;
            my @f = split /\t/, $line, -1;
            next unless @f >= 2;
            my $bp = 0 + $f[0];
            my $y = 0 + $f[1];
            my $x = locus_bp_to_column_x($i, $start, $end, $bp);
            print {$out} join("\t", sprintf('%.6f', $x), sprintf('%.4f', $y), $i + 1), "\n";
        }
        close $fh;
    }
    close $out or die "Cannot close $combined_tsv: $!\n";

    my $top_logp = $args{top_logp} || 8;
    my $track_count = scalar @{ $args{labels} || [] };
    my $ymax = $track_count * $top_logp + 1.2;
    my $ymin = ($annotation_mode eq 'gtf') ? -8.8 : -10.5;
    my $xmax = $n;
    my $title = $args{title} || 'Local top hits Manhattan Plot';
    my $sig_y = $args{sig_y};
    $sig_y = safe_neglog10_text('1e-6') unless defined $sig_y;

    open my $gp, '>', $gp_file or die "Cannot write $gp_file: $!\n";
    print {$gp} "set terminal png enhanced size 2200,1780\n";
    print {$gp} "set output '" . escape_gp($args{output_png}) . "'\n";
    print {$gp} "set datafile separator '\\t'\n";
    print {$gp} "set title \"" . escape_gp($title) . "\"\n";
    print {$gp} "set xrange [0:$xmax]\n";
    print {$gp} "set yrange [$ymin:$ymax]\n";
    print {$gp} "set xlabel ''\n";
    print {$gp} "set ylabel '-Log10(p)'\n";
    print {$gp} "set border 3\n";
    print {$gp} "set lmargin 8\n";
    print {$gp} "set rmargin 3\n";
    print {$gp} "set bmargin 13\n";
    print {$gp} "unset key\n";
    print {$gp} "set tics out nomirror\n";
    print {$gp} "set grid ytics lc rgb '#dddddd' dt 2\n";
    print {$gp} "unset xtics\n";
    print {$gp} "set ytics (" . join(', ', repeated_panel_ytics_wrapper($track_count, $top_logp)) . ")\n";
    print {$gp} "set mytics 2\n";

    my $arrow_id = 10;
    for my $track_i (0 .. ($track_count - 1)) {
        my $base = $track_i * $top_logp;
        print {$gp} "set arrow $arrow_id from 0," . ($base + $sig_y) . " to $xmax," . ($base + $sig_y) . " nohead dt 2 lc rgb '#777777' lw 1\n";
        $arrow_id++;
        next unless $track_i > 0;
        print {$gp} "set arrow $arrow_id from 0,$base to $xmax,$base nohead lc rgb '#bbbbbb' lw 1\n";
        $arrow_id++;
    }
    for my $i (1 .. ($n - 1)) {
        print {$gp} "set arrow $arrow_id from $i,$ymin to $i,$ymax nohead lc rgb '#c8c8c8' lw 1\n";
        $arrow_id++;
    }
    for my $track_i (0 .. ($track_count - 1)) {
        my $panel_y = $track_i * $top_logp + $top_logp - 0.95;
        print {$gp} "set label " . (100 + $track_i) . " \"" . escape_gp($args{labels}[$track_i]) . "\" at " . ($xmax / 2) . ",$panel_y center font ',18'\n";
    }

    my $label_id = 1000;
    if ($annotation_mode eq 'gtf') {
        my $obj_id = 4000;
        for my $i (0 .. $#items) {
            my $center = $i + 0.5;
            print {$gp} "set arrow $arrow_id from $center,0 to $center,$ymax nohead lc rgb '#a0a0a0' dt 2 lw 1\n";
            $arrow_id++;
            my ($genes, $lane_count) = mini_gtf_genes_for_column($items[$i]);
            next unless @{$genes};
            my $track_top = -0.95;
            my $lane_gap = 1.20;
            my $bar_h = 0.22;
            for my $g (@{$genes}) {
                my $lane = $g->{lane} || 0;
                my $bar_top = $track_top - $lane * $lane_gap;
                my $bar_bot = $bar_top - $bar_h;
                my $label_y = $bar_bot - 0.16;
                print {$gp} "set object $obj_id rect from $g->{x1},$bar_bot to $g->{x2},$bar_top fc rgb '$g->{color}' fillstyle solid 0.72 border lc rgb '$g->{color}' lw 1 front\n";
                $obj_id++;
                print {$gp} "set label $label_id \"" . escape_gp($g->{gene}) . "\" at $g->{xc},$label_y center font '" . italic_font_spec_gp(10) . "'\n";
                $label_id++;
            }
        }
    }
    elsif ($annotation_mode eq 'labels') {
        for my $i (0 .. $#items) {
            my $center = $i + 0.5;
            my $label_dx = 0.045;
            my $label_y = $ymin * 0.48;
            my $snp = $items[$i]{snp} || '';
            my $gene = nearest_gene_label_for_item($items[$i]);
            if (!length $gene) {
                my $manifest = read_manifest_tsv($items[$i]{manifest});
                my $bp = $items[$i]{bp} || $manifest->{TARGET_BP} || $manifest->{target_bp};
                my $chr = $items[$i]{chr} || $manifest->{TARGET_CHR} || $manifest->{target_chr};
                my $gtf_file = $items[$i]{gtf_file} || infer_cached_locus_gtf_path(
                    output_dir => $args{output_dir},
                    snp        => $items[$i]{snp},
                    window_bp  => ($manifest->{WINDOW_BP} || $manifest->{window_bp} || 0),
                    runner     => {},
                );
                $gene = nearest_gene_from_gtf_file(
                    gtf_file => $gtf_file,
                    chr      => $chr,
                    bp       => $bp,
                ) if $gtf_file;
                if (!length $gene && $gtf_file) {
                    $gene = nearest_gene_from_gtf_file(
                        gtf_file => $gtf_file,
                        chr      => '',
                        bp       => $bp,
                    );
                }
                if (!length $gene) {
                    $gene = nearest_gene_for_snp_from_cached_gtf(
                        output_dir => $args{output_dir},
                        snp        => $items[$i]{snp},
                        bp         => $bp,
                        window_bp  => ($manifest->{WINDOW_BP} || $manifest->{window_bp} || 0),
                    );
                }
            }
            print {$gp} "set arrow $arrow_id from $center,0 to $center,$ymax nohead lc rgb '#a0a0a0' dt 2 lw 1\n";
            $arrow_id++;
            if (length $gene) {
                print {$gp} "set label $label_id \"" . escape_gp($snp) . "\" at " . ($center - $label_dx) . ",$label_y center rotate by 90 font ',12'\n";
                $label_id++;
                print {$gp} "set label $label_id \"" . escape_gp($gene) . "\" at " . ($center + $label_dx) . ",$label_y center rotate by 90 font '" . italic_font_spec_gp(12) . "'\n";
                $label_id++;
            }
            else {
                print {$gp} "set label $label_id \"" . escape_gp($snp) . "\" at $center,$label_y center rotate by 90 font ',12'\n";
                $label_id++;
            }
        }
    }

    print {$gp} "plot ";
    my @plots;
    for my $i (0 .. $#items) {
        my $color = chromosome_color_wrapper($items[$i]{chr});
        push @plots,
            "'" . escape_gp($combined_tsv) . "' using ((\$3==" . ($i + 1) . ")?\$1:1/0):2 with points pt 7 ps 0.9 lc rgb '$color'";
    }
    print {$gp} join(", \\\n     ", @plots) . "\n";
    close $gp or die "Cannot close $gp_file: $!\n";

    system($args{gnuplot}, $gp_file) == 0
        or die "gnuplot failed for combined local Manhattan batch $gp_file\n";
}

sub render_combined_local_manhattan_gtf_batch {
    my (%args) = @_;
    my @items = @{ $args{images} || [] };
    die "No local Manhattan loci supplied for combined GTF batch render\n" unless @items;
    my $n = scalar @items;
    my $cols = $args{columns} || $n;
    $cols = $n if $cols > $n;
    $cols = 1 if $cols < 1;

    my $combined_scaled_tsv = File::Spec->catfile($args{output_dir}, $args{output_base} . '.combined_scaled.tsv');
    my $gp_file = File::Spec->catfile($args{output_dir}, $args{output_base} . '.combined_gtf.gp');

    open my $scaled_out, '>', $combined_scaled_tsv or die "Cannot write $combined_scaled_tsv: $!\n";
    print {$scaled_out} join("\t", qw(KIND KIND_CODE LOCUS X1 X2 Y TRACK COLORVAL GENE LANE BP SNP IS_TARGET)), "\n";

    my $max_lane = 0;
    for my $i (0 .. $#items) {
        my $item = $items[$i];
        my $manifest = read_manifest_tsv($item->{manifest});
        my $target_bp = $manifest->{BP} || $manifest->{bp} || $manifest->{TARGET_BP} || $manifest->{target_bp} || $item->{bp};
        my $window_bp = $manifest->{WINDOW_BP} || $manifest->{window_bp} || 0;
        my $start = $manifest->{WINDOW_START} || $manifest->{window_start} || ($target_bp - $window_bp);
        my $end   = $manifest->{WINDOW_END} || $manifest->{window_end} || ($target_bp + $window_bp);
        $start = 1 if !defined $start || $start < 1;
        $end = $start + 1 if !defined $end || $end <= $start;

        open my $pfh, '<', $item->{plot_tsv} or die "Cannot read $item->{plot_tsv}: $!\n";
        my $pheader = <$pfh>;
        while (my $line = <$pfh>) {
            chomp $line;
            $line =~ s/\r$//;
            next unless length $line;
            my @f = split /\t/, $line, -1;
            next unless @f >= 7;
            my $bp = 0 + $f[0];
            my $x = locus_bp_to_column_x($i, $start, $end, $bp);
            print {$scaled_out} join(
                "\t",
                'assoc',
                1,
                $i + 1,
                sprintf('%.6f', $x),
                sprintf('%.6f', $x),
                $f[1],
                $f[2],
                $f[6],
                '',
                '',
                $f[0],
                $f[5],
                $f[4],
            ), "\n";
        }
        close $pfh;

        next unless $item->{manifest} && -f $item->{manifest};
        my $gene_tsv = $manifest->{gene_tsv} || $manifest->{GENE_TSV} || '';
        $gene_tsv = $item->{image} if !$gene_tsv;
        $gene_tsv =~ s/\.png$/.genes.tsv/i if $gene_tsv;
        next unless $gene_tsv && -f $gene_tsv;
        open my $gfh, '<', $gene_tsv or die "Cannot read $gene_tsv: $!\n";
        my $gheader = <$gfh>;
        while (my $line = <$gfh>) {
            chomp $line;
            $line =~ s/\r$//;
            next unless length $line;
            my ($type, $gene, $gstart, $gend, $lane) = split /\t/, $line, -1;
            next unless defined $type && $type =~ /^(?:gene|exon)$/;
            next unless defined $gstart && defined $gend && $gstart =~ /^\d+$/ && $gend =~ /^\d+$/;
            my $x1 = locus_bp_to_column_x($i, $start, $end, $gstart);
            my $x2 = locus_bp_to_column_x($i, $start, $end, $gend);
            my ($col_min, $col_max) = locus_column_bounds($i);
            $x1 = $col_min if $x1 < $col_min;
            $x2 = $col_max if $x2 > $col_max;
            if ($x2 <= $x1) {
                my $min_w = 0.010;
                my $mid = ($x1 + $x2) / 2;
                $x1 = $mid - $min_w / 2;
                $x2 = $mid + $min_w / 2;
                $x1 = $col_min if $x1 < $col_min;
                $x2 = $col_max if $x2 > $col_max;
            }
            print {$scaled_out} join(
                "\t",
                $type,
                ($type eq 'gene' ? 2 : 3),
                $i + 1,
                sprintf('%.6f', $x1),
                sprintf('%.6f', $x2),
                '',
                '',
                '',
                $gene,
                $lane,
                $gstart,
                '',
                '',
            ), "\n";
            $max_lane = $lane if defined $lane && $lane =~ /^\d+$/ && $lane > $max_lane;
        }
        close $gfh;
    }
    close $scaled_out or die "Cannot close $combined_scaled_tsv: $!\n";

    my $top_logp = $args{top_logp} || 8;
    my $track_count = scalar @{ $args{labels} || [] };
    my $ymax = $track_count * $top_logp + 1.2;
    my $gene_height = 1.8 * ($max_lane + 1);
    $gene_height = 6.0 if $gene_height < 6.0;
    my $ymin = -$gene_height;
    my $xmax = $n;
    my $title = $args{title} || 'Local top hits Manhattan Plot';
    my $sig_y = $args{sig_y};
    $sig_y = safe_neglog10_text('1e-6') unless defined $sig_y;

    open my $gp, '>', $gp_file or die "Cannot write $gp_file: $!\n";
    print {$gp} "set terminal png enhanced size 2200,1780\n";
    print {$gp} "set output '" . escape_gp($args{output_png}) . "'\n";
    print {$gp} "set datafile separator '\\t'\n";
    print {$gp} "set title \"" . escape_gp($title) . "\"\n";
    print {$gp} "set xrange [0:$xmax]\n";
    print {$gp} "set yrange [$ymin:$ymax]\n";
    print {$gp} "set xlabel ''\n";
    print {$gp} "set ylabel '-Log10(p)'\n";
    print {$gp} "set border 3\n";
    print {$gp} "set lmargin 8\n";
    print {$gp} "set rmargin 3\n";
    print {$gp} "set bmargin 13\n";
    print {$gp} "unset key\n";
    print {$gp} "set tics out nomirror\n";
    print {$gp} "set grid ytics lc rgb '#dddddd' dt 2\n";
    print {$gp} "unset xtics\n";
    print {$gp} "set ytics (" . join(', ', repeated_panel_ytics_wrapper($track_count, $top_logp)) . ")\n";
    print {$gp} "set mytics 2\n";

    my $arrow_id = 10;
    for my $track_i (0 .. ($track_count - 1)) {
        my $base = $track_i * $top_logp;
        print {$gp} "set arrow $arrow_id from 0," . ($base + $sig_y) . " to $xmax," . ($base + $sig_y) . " nohead dt 2 lc rgb '#777777' lw 1\n";
        $arrow_id++;
        next unless $track_i > 0;
        print {$gp} "set arrow $arrow_id from 0,$base to $xmax,$base nohead lc rgb '#bbbbbb' lw 1\n";
        $arrow_id++;
    }
    for my $i (1 .. ($n - 1)) {
        print {$gp} "set arrow $arrow_id from $i,$ymin to $i,$ymax nohead lc rgb '#c8c8c8' lw 1\n";
        $arrow_id++;
    }
    for my $track_i (0 .. ($track_count - 1)) {
        my $panel_y = $track_i * $top_logp + $top_logp - 0.95;
        print {$gp} "set label " . (100 + $track_i) . " \"" . escape_gp($args{labels}[$track_i]) . "\" at " . ($xmax / 2) . ",$panel_y center font ',18'\n";
    }

    my @gene_colors = (
        '#f4a6a6', '#f7c97f', '#d8d35f', '#9fd27a', '#71c9b8',
        '#86b6f6', '#a995e8', '#de9de6', '#f4a8c4', '#c8b08b',
        '#ef8b62', '#9ec3a5', '#7fb3d5', '#c39bd3', '#f8c471',
    );
    open my $gfh, '<', $combined_scaled_tsv or die "Cannot read $combined_scaled_tsv: $!\n";
    <$gfh>;
    my %gene_color_for;
    my %lane_label_count;
    my $obj_id = 4000;
    my $label_id = 5000;
    while (my $line = <$gfh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my ($type, $kind_code, $locus, $x1, $x2, $y, $track, $colorval, $gene, $lane) = split /\t/, $line, -1;
        next unless defined $kind_code && $kind_code =~ /^(?:2|3)$/;
        next unless defined $type && defined $gene;
        my $track_top = -0.55 - (1.35 * $lane);
        my $line_top = $track_top - 0.12;
        my $line_bot = $line_top - 0.22;
        my $exon_top = $track_top + 0.02;
        my $exon_bot = $exon_top - 0.52;
        my $mid = ($x1 + $x2) / 2;
        my $gene_key = join('|', $locus, $gene);
        my $gene_color = $gene_color_for{$gene_key};
        if (!defined $gene_color) {
            my $color_idx = scalar(keys %gene_color_for) % @gene_colors;
            $gene_color = $gene_colors[$color_idx];
            $gene_color_for{$gene_key} = $gene_color;
        }
        if ($type eq 'gene') {
            my $lane_key = join('|', $locus, $lane);
            my $lane_rank = $lane_label_count{$lane_key} // 0;
            my @below_label_offsets = (0.08, 0.22, 0.36, 0.50);
            my $label_y = $line_bot - $below_label_offsets[$lane_rank % @below_label_offsets];
            $lane_label_count{$lane_key} = $lane_rank + 1;
            print {$gp} "set object $obj_id rect from $x1,$line_bot to $x2,$line_top fc rgb '$gene_color' fillstyle solid 0.72 border lc rgb '$gene_color' lw 1 behind\n";
            $obj_id++;
            print {$gp} "set label $label_id \"" . escape_gp($gene) . "\" at $mid,$label_y center front tc rgb '#111111' font '" . italic_font_spec_gp(10) . "'\n";
            $label_id++;
        } else {
            my $disp_start = $x1;
            my $disp_end = $x2;
            my $min_w = 0.010;
            if (($disp_end - $disp_start) < $min_w) {
                my $half_w = $min_w / 2;
                $disp_start = $mid - $half_w;
                $disp_end = $mid + $half_w;
            }
            print {$gp} "set object $obj_id rect from $disp_start,$exon_bot to $disp_end,$exon_top fc rgb '$gene_color' fillstyle solid 1.0 border lc rgb '$gene_color' lw 1 front\n";
            $obj_id++;
        }
    }
    close $gfh;

    my @plots = (
        "'" . escape_gp($combined_scaled_tsv) . "' using ((\$2==1)?\$4:1/0):6:8 with points pt 7 ps 0.72 lc palette"
    );
    print {$gp} "set palette maxcolors 12 defined (1 '#1f77b4', 2 '#ff7f0e', 3 '#2ca02c', 4 '#d62728', 5 '#9467bd', 6 '#8c564b', 7 '#e377c2', 8 '#7f7f7f', 9 '#bcbd22', 10 '#17becf', 11 '#3366cc', 12 '#dd4477')\n";
    print {$gp} "unset colorbox\n";
    print {$gp} "plot " . join(", \\\n     ", @plots) . "\n";
    close $gp or die "Cannot close $gp_file: $!\n";

    system($args{gnuplot}, $gp_file) == 0
        or die "gnuplot failed for combined local Manhattan GTF batch $gp_file\n";
}

sub prepare_locus_wide_sources {
    my (%args) = @_;
    my %sources;
    my $safe_window = safe_name($args{window_bp});

    for my $hit (@{ $args{hits} || [] }) {
        my $snp = $hit->{SNP} || '';
        next unless length $snp;
        my $safe_snp = safe_name($snp);
        if (!defined($hit->{CHR}) || !defined($hit->{BP})) {
            for my $coordinate_manifest (glob(File::Spec->catfile($args{output_dir}, '*' . $safe_snp . '*.manifest.tsv'))) {
                next unless -s $coordinate_manifest;
                my $coordinate_metrics = read_manifest_tsv($coordinate_manifest);
                my $coordinate_snp = $coordinate_metrics->{target_snp} // $coordinate_metrics->{snp} // '';
                next unless lc($coordinate_snp) eq lc($snp);
                my $coordinate_chr = $coordinate_metrics->{target_chr} // $coordinate_metrics->{chr};
                my $coordinate_bp  = $coordinate_metrics->{target_bp}  // $coordinate_metrics->{bp};
                next unless defined $coordinate_chr && $coordinate_chr ne '' && defined $coordinate_bp && $coordinate_bp ne '';
                $hit->{CHR} = $coordinate_chr if !defined $hit->{CHR};
                $hit->{BP}  = $coordinate_bp  if !defined $hit->{BP};
                last;
            }
        }
        my $data = File::Spec->catfile(
            $args{output_dir},
            "gunplot_locus_${safe_snp}_window_${safe_window}.wide.tsv.gz",
        );
        my $manifest = File::Spec->catfile(
            $args{output_dir},
            "gunplot_locus_${safe_snp}_window_${safe_window}.wide.manifest.tsv",
        );
        my ($valid, $metrics) = locus_wide_cache_matches(
            data       => $data,
            manifest   => $manifest,
            snp        => $snp,
            window_bp  => $args{window_bp},
            exact_window => 1,
        );

        if (($args{force} || !$valid)
            && $args{source_long} && -s $args{source_long}
            && $args{preset_config} && -f $args{preset_config}) {
            my @extract_cmd = (
                $^X,
                File::Spec->catfile($Bin, 'DiffGWASDeps', 'extract_single_snp_wide_diff_gwas.pl'),
                '--config', $args{preset_config},
                '--input', $args{source_long},
                '--target-snp', $snp,
                '--window-bp', $args{window_bp},
                '--output', $data,
                '--manifest', $manifest,
                '--output-dir', $args{output_dir},
            );
            push @extract_cmd, ('--target-chr', $hit->{CHR}) if defined $hit->{CHR} && $hit->{CHR} ne '';
            push @extract_cmd, ('--target-bp',  $hit->{BP})  if defined $hit->{BP}  && $hit->{BP}  ne '';
            run_cmd(\@extract_cmd, "tabix-backed metadata locus for $snp");
            ($valid, $metrics) = locus_wide_cache_matches(
                data       => $data,
                manifest   => $manifest,
                snp        => $snp,
                window_bp  => $args{window_bp},
                exact_window => 1,
            );
            die "Target-locus extraction did not produce a valid cache for $snp\n" unless $valid;
        }

        if (!$valid) {
            my @compatible;
            for my $candidate (glob(File::Spec->catfile($args{output_dir}, '*' . $safe_snp . '.wide.tsv.gz'))) {
                (my $candidate_manifest = $candidate) =~ s/\.wide\.tsv\.gz$/.wide.manifest.tsv/;
                my ($candidate_valid, $candidate_metrics) = locus_wide_cache_matches(
                    data       => $candidate,
                    manifest   => $candidate_manifest,
                    snp        => $snp,
                    window_bp  => $args{window_bp},
                    exact_window => 0,
                );
                next unless $candidate_valid;
                push @compatible, [$candidate, $candidate_manifest, $candidate_metrics];
            }
            if (@compatible) {
                @compatible = sort {
                    (0 + ($a->[2]{window_bp} || 9e99)) <=> (0 + ($b->[2]{window_bp} || 9e99))
                } @compatible;
                ($data, $manifest, $metrics) = @{ $compatible[0] };
                $valid = 1;
            }
        }

        next unless $valid;
        $hit->{CHR} = $metrics->{target_chr} if !defined($hit->{CHR}) && defined($metrics->{target_chr});
        $hit->{BP}  = $metrics->{target_bp}  if !defined($hit->{BP})  && defined($metrics->{target_bp});
        $sources{uc($snp)} = {
            data     => $data,
            manifest => $manifest,
        };
        print "[prep] Using compact target-locus metadata source for $snp: $data\n";
    }
    return \%sources;
}

sub locus_wide_cache_matches {
    my (%args) = @_;
    return (0, {}) unless $args{data} && -s $args{data} && $args{manifest} && -s $args{manifest};
    my $metrics = read_manifest_tsv($args{manifest});
    my $manifest_snp = $metrics->{target_snp} // $metrics->{snp} // '';
    return (0, $metrics) unless lc($manifest_snp) eq lc($args{snp} || '');
    my $manifest_window = $metrics->{window_bp};
    return (0, $metrics) unless defined $manifest_window && $manifest_window ne '';
    my $window_ok = $args{exact_window}
        ? ((0 + $manifest_window) == (0 + $args{window_bp}))
        : ((0 + $manifest_window) >= (0 + $args{window_bp}));
    return ($window_ok ? 1 : 0, $metrics);
}

sub write_top_hits_csv {
    my (%args) = @_;
    my ($existing_header, $existing_rows) = extract_existing_sas_style_top_hits_csv_rows(
        runner => $args{runner},
        hits   => $args{hits},
    );
    my ($fallback_header, $fallback_rows) = (!$existing_header || !$existing_rows)
        ? build_sas_style_top_hits_from_wide(
            runner    => $args{runner},
            wide_data => $args{wide_data},
            hits      => $args{hits},
            locus_sources => $args{locus_sources},
          )
        : (undef, undef);
    my @header = $existing_header ? @{$existing_header} : @{$fallback_header || []};

    open my $out, '>', $args{path} or die "Cannot write $args{path}: $!\n";
    print {$out} join(',', map { csv_escape($_) } @header), "\n";
    for my $i (0 .. $#{ $args{hits} }) {
        my $hit = $args{hits}[$i];
        my $panel_index = int($i / ($args{batch_size} || 1)) + 1;
        my $row = $existing_rows ? $existing_rows->{ uc($hit->{SNP} || '') } : $fallback_rows->{ uc($hit->{SNP} || '') };
        my %values = $row ? %{$row} : ();
        $values{CHR} = $hit->{CHR} if defined $hit->{CHR} && length $hit->{CHR};
        $values{BP} = $hit->{BP} if defined $hit->{BP} && length $hit->{BP};
        $values{SNP} = $hit->{SNP} if defined $hit->{SNP} && length $hit->{SNP};
        if (defined $hit->{gene} && length $hit->{gene}) {
            $values{gene} = $hit->{gene};
            $values{snp_gene} = $hit->{snp_gene} if defined $hit->{snp_gene} && length $hit->{snp_gene};
            $values{gene_source} = defined $hit->{gene_source} && length $hit->{gene_source}
                ? $hit->{gene_source}
                : 'USER';
        }
        $values{hit_order} = $i + 1;
        $values{panel_index} = $panel_index;
        print {$out} join(',', map { csv_escape(defined $values{$_} ? $values{$_} : '') } @header), "\n";
    }
    close $out or die "Cannot close $args{path}: $!\n";
}

sub extract_existing_sas_style_top_hits_csv_rows {
    my (%args) = @_;
    my $runner = $args{runner} || {};
    my %wanted = map { uc($_->{SNP} || '') => 1 } @{ $args{hits} || [] };
    my $csv_name = $runner->{LOCAL_TOP_HITS_CSV_BASENAME} || '';
    return unless $csv_name;
    my @candidate_paths = (
        File::Spec->catfile(localize_path($runner->{OUTPUT_DIR} || ''), $csv_name),
        File::Spec->catfile($Bin, $csv_name),
        $csv_name,
    );
    my ($path) = grep { defined $_ && length $_ && -f $_ } @candidate_paths;
    return unless $path;
    open my $fh, '<', $path or return;
    my $header = <$fh>;
    return unless defined $header;
    chomp $header;
    $header =~ s/\r$//;
    my @cols = map { trim($_) } parse_csv_line($header);
    return unless @cols;
    my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
    return unless exists $idx{SNP};
    my %rows_by_snp;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = parse_csv_line($line);
        my %row;
        @row{@cols} = @f;
        my $snp = uc($row{SNP} || '');
        next unless $wanted{$snp};
        $rows_by_snp{$snp} = \%row;
    }
    close $fh;
    return (\@cols, \%rows_by_snp);
}

sub build_sas_style_top_hits_from_wide {
    my (%args) = @_;
    my $runner = $args{runner} || {};
    my $wide_data = $args{wide_data} || '';
    return ([], {}) unless $wide_data && -s $wide_data;

    my %targets = map {
        my $snp = uc($_->{SNP} || '');
        $snp => $_
    } @{ $args{hits} || [] };
    my %rows_by_snp;
    my @wide_cols;
    my %idx;
    my $locus_sources = $args{locus_sources} || {};
    for my $snp (sort keys %targets) {
        my $source = $locus_sources->{$snp};
        next unless $source && $source->{data} && -s $source->{data};
        my ($cols, $row) = read_target_row_from_wide($source->{data}, $snp);
        next unless $cols && $row;
        @wide_cols = @$cols unless @wide_cols;
        $rows_by_snp{$snp} = $row;
    }

    my %missing = map { $_ => 1 } grep { !exists $rows_by_snp{$_} } keys %targets;
    if (%missing) {
        print "[prep] Compact locus metadata unavailable for " . scalar(keys %missing)
            . " target(s); scanning the full wide table as fallback.\n";
        my $fh = IO::Uncompress::Gunzip->new($wide_data)
            or die "Cannot read $wide_data: $GunzipError\n";
        my $header = <$fh>;
        if (defined $header) {
            chomp $header;
            $header =~ s/\r$//;
            @wide_cols = split /\t/, $header, -1;
            %idx = map { $wide_cols[$_] => $_ } 0 .. $#wide_cols;
        }
        while (my $line = <$fh>) {
            chomp $line;
            $line =~ s/\r$//;
            next unless length $line;
            my @f = split /\t/, $line, -1;
            my $snp = uc($f[$idx{SNP}] // '');
            next unless $missing{$snp};
            my %row;
            @row{@wide_cols} = @f;
            $rows_by_snp{$snp} = \%row;
            delete $missing{$snp};
            last unless %missing;
        }
        close $fh;
    }
    else {
        print "[prep] Top-hit CSV metadata loaded from compact target-locus data; full wide-table scan skipped.\n";
    }

    my @header = (
        qw(hit_order panel_index CHR BP SNP EFFECT_ALLELE OTHER_ALLELE REFERENCE_ALLELE ALTERNATIVE_ALLELE gene snp_gene focus_signal selected_maf maf_source gwas_group1_maf gwas_group2_maf gwas_pair_maf_min gnomad_maf gnomad_pops maf_filter_decision maf_filter_reason)
    );
    my %skip = map { $_ => 1 } qw(CHR BP SNP gene snp_gene A1 A2 REQUESTED_HIT_ORDER hit_order panel_index focus_signal EFFECT_ALLELE OTHER_ALLELE REFERENCE_ALLELE ALTERNATIVE_ALLELE selected_maf maf_source gwas_group1_maf gwas_group2_maf gwas_pair_maf_min gnomad_maf gnomad_pops maf_filter_decision maf_filter_reason);
    my @extra_cols = grep { !$skip{$_} } @wide_cols;
    push @header, @extra_cols;
    push @header, 'gene_source' unless grep { $_ eq 'gene_source' } @header;

    my %export_rows;
    for my $hit (@{ $args{hits} || [] }) {
        my $snp = uc($hit->{SNP} || '');
        my %wide = %{ $rows_by_snp{$snp} || {} };
        my %row = %wide;
        my $gene = $hit->{gene};
        $gene = extract_gene_from_snp_gene($hit->{snp_gene}) if !defined($gene) || !length($gene);
        $gene = '' unless defined $gene;
        my $snp_gene = $hit->{snp_gene};
        $snp_gene = ($hit->{SNP} || '') . ':' . (length($gene) ? $gene : 'NA')
            unless defined($snp_gene) && length($snp_gene);
        my $focus_pvar = $runner->{TOP_HIT_FOCUS_PVAR} || '';
        my $focus_signal = (length($focus_pvar) && exists $wide{$focus_pvar}) ? $wide{$focus_pvar} : $hit->{focus_signal};
        $row{CHR} = defined $hit->{CHR} && length($hit->{CHR}) ? $hit->{CHR} : $wide{CHR};
        $row{BP}  = defined $hit->{BP}  && length($hit->{BP})  ? $hit->{BP}  : $wide{BP};
        $row{SNP} = $hit->{SNP};
        $row{EFFECT_ALLELE} = defined $wide{A1} ? $wide{A1} : '';
        $row{OTHER_ALLELE} = defined $wide{A2} ? $wide{A2} : '';
        $row{REFERENCE_ALLELE} = defined $wide{A2} ? $wide{A2} : '';
        $row{ALTERNATIVE_ALLELE} = defined $wide{A1} ? $wide{A1} : '';
        $row{gene} = $gene;
        $row{snp_gene} = $snp_gene;
        $row{focus_signal} = defined $focus_signal ? $focus_signal : '';
        $row{selected_maf} = $hit->{selected_maf} if defined $hit->{selected_maf} && length $hit->{selected_maf};
        $row{maf_source} = $hit->{maf_source} if defined $hit->{maf_source} && length $hit->{maf_source};
        $row{gwas_group1_maf} = $hit->{gwas_group1_maf} if defined $hit->{gwas_group1_maf} && length $hit->{gwas_group1_maf};
        $row{gwas_group2_maf} = $hit->{gwas_group2_maf} if defined $hit->{gwas_group2_maf} && length $hit->{gwas_group2_maf};
        $row{gwas_pair_maf_min} = $hit->{gwas_pair_maf_min} if defined $hit->{gwas_pair_maf_min} && length $hit->{gwas_pair_maf_min};
        $row{gnomad_maf} = $hit->{gnomad_maf} if defined $hit->{gnomad_maf} && length $hit->{gnomad_maf};
        $row{gnomad_pops} = $hit->{gnomad_pops} if defined $hit->{gnomad_pops} && length $hit->{gnomad_pops};
        $row{maf_filter_decision} = $hit->{maf_filter_decision} if defined $hit->{maf_filter_decision} && length $hit->{maf_filter_decision};
        $row{maf_filter_reason} = $hit->{maf_filter_reason} if defined $hit->{maf_filter_reason} && length $hit->{maf_filter_reason};
        $row{gene_source} = defined $hit->{gene_source} && length($hit->{gene_source})
            ? $hit->{gene_source}
            : (defined $wide{gene_source} ? $wide{gene_source} : (length($gene) ? 'NA' : 'NA'));
        $export_rows{$snp} = \%row;
    }
    return (\@header, \%export_rows);
}

sub read_target_row_from_wide {
    my ($path, $target_snp) = @_;
    return unless $path && -s $path && defined $target_snp && length $target_snp;
    my $fh = IO::Uncompress::Gunzip->new($path)
        or die "Cannot read $path: $GunzipError\n";
    my $header = <$fh>;
    return unless defined $header;
    chomp $header;
    $header =~ s/\r$//;
    my @cols = split /\t/, $header, -1;
    my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
    unless (exists $idx{SNP}) {
        close $fh;
        return;
    }
    my $wanted = uc($target_snp);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        next unless uc($f[$idx{SNP}] // '') eq $wanted;
        my %row;
        @row{@cols} = @f;
        close $fh;
        return (\@cols, \%row);
    }
    close $fh;
    return;
}

sub parse_csv_line {
    my ($line) = @_;
    my @fields = parse_line(',', 0, $line);
    return @fields;
}

sub read_forest_manifest_rows {
    my ($path) = @_;
    return () unless defined $path && length $path && -f $path;
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    my $header = <$fh>;
    return () unless defined $header;
    chomp $header;
    $header =~ s/\r$//;
    my @cols = split /\t/, $header, -1;
    my @rows;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        my %row;
        @row{@cols} = @f;
        push @rows, \%row;
    }
    close $fh;
    return @rows;
}

sub forest_csv_basename_for_targets {
    my (%args) = @_;
    my $name = $args{base_name} || 'gunplot_forest_top_hits.csv';
    my $targets = $args{target_snps} || '';
    return $name unless length $targets;
    my @items = grep { length } map { trim($_) } split /,/, $targets;
    return $name unless @items;
    my $variant = @items == 1 ? safe_name($items[0]) : ('targets_' . scalar(@items));
    if ($name =~ /^(.*?)(\.[^.]+)$/) {
        return $1 . '_' . $variant . $2;
    }
    return $name . '_' . $variant;
}

sub write_single_image_html {
    my (%args) = @_;
    open my $out, '>', $args{path} or die "Cannot write $args{path}: $!\n";
    print {$out} <<"HTML";
<!doctype html>
<html><head><meta charset="utf-8">
<title>@{[ html_escape($args{title}) ]}</title>
<style>body{margin:0;padding:16px;font-family:Arial,sans-serif;background:#fff;} img{max-width:100%;height:auto;display:block;}</style>
</head><body>
<img src="@{[ html_escape($args{image_rel}) ]}" alt="@{[ html_escape($args{title}) ]}">
</body></html>
HTML
    close $out or die "Cannot close $args{path}: $!\n";
}

sub write_gallery_html {
    my (%args) = @_;
    open my $out, '>', $args{path} or die "Cannot write $args{path}: $!\n";
    print {$out} <<"HTML";
<!doctype html>
<html><head><meta charset="utf-8"><title>@{[ html_escape($args{title}) ]}</title>
<style>
body{font-family:Arial,Helvetica,sans-serif;margin:24px;background:#fff;}
img{max-width:100%;height:auto;display:block;border:1px solid #ddd;margin:8px 0 24px 0;}
h1{font-size:20px}
h2{font-size:15px;margin-top:20px}
a{color:#0b61a4}
</style></head><body>
<h1>@{[ html_escape($args{title}) ]}</h1>
HTML
    if ($args{top_csv}) {
        print {$out} '<p><a href="' . html_escape($args{top_csv}) . '">Top-hit CSV</a></p>' . "\n";
    }
    for my $img (@{ $args{images} }) {
        my $caption = $img->{snp} || 'locus';
        if (defined $img->{chr} && $img->{chr} ne '' && defined $img->{bp} && $img->{bp} ne '') {
            $caption .= " ($img->{chr}:$img->{bp})";
        }
        print {$out} '<h2>' . html_escape($caption) . "</h2>\n";
        print {$out} '<img src="' . html_escape(basename($img->{image})) . '" alt="' . html_escape($caption) . qq{">\n};
    }
    print {$out} "</body></html>\n";
    close $out or die "Cannot close $args{path}: $!\n";
}

sub write_index_html {
    my (%args) = @_;
    open my $out, '>', $args{path} or die "Cannot write $args{path}: $!\n";
    print {$out} <<"HTML";
<!doctype html>
<html><head><meta charset="utf-8"><title>@{[ html_escape($args{title}) ]}</title></head>
<body style="font-family:Arial,Helvetica,sans-serif;margin:24px">
<h1 style="font-size:20px">@{[ html_escape($args{title}) ]}</h1>
<p>@{[ html_escape($args{batch_msg}) ]}</p>
<ul>
HTML
    for my $part (@{ $args{parts} }) {
        print {$out} '<li><a href="' . html_escape(basename($part)) . '">' . html_escape(basename($part)) . "</a></li>\n";
    }
    print {$out} "</ul>\n";
    if ($args{top_csv}) {
        print {$out} '<p><a href="' . html_escape($args{top_csv}) . '">Combined top-hit CSV</a></p>' . "\n";
    }
    print {$out} "</body></html>\n";
    close $out or die "Cannot close $args{path}: $!\n";
}

sub read_hits_tsv {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    my $header = <$fh>;
    die "Top-hit table is empty: $path\n" unless defined $header;
    chomp $header;
    my @cols = split /\t/, $header, -1;
    my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
    my @hits;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        push @hits, {
            hit_order    => pick_field(\%idx, \@f, qw(hit_order locus_rank candidate_rank)),
            CHR          => pick_field(\%idx, \@f, qw(CHR chr)),
            BP           => pick_field(\%idx, \@f, qw(BP bp)),
            SNP          => pick_field(\%idx, \@f, qw(SNP snp rsid)),
            gene         => pick_field(\%idx, \@f, qw(gene nearest_gene genesymbol)),
            snp_gene     => pick_field(\%idx, \@f, qw(snp_gene SNP_GENE)),
            focus_signal => pick_field(\%idx, \@f, qw(focus_signal common_assoc_p)),
            selected_maf => pick_field(\%idx, \@f, qw(selected_maf)),
            maf_source   => pick_field(\%idx, \@f, qw(maf_source)),
            gwas_group1_maf => pick_field(\%idx, \@f, qw(gwas_group1_maf)),
            gwas_group2_maf => pick_field(\%idx, \@f, qw(gwas_group2_maf)),
            gwas_pair_maf_min => pick_field(\%idx, \@f, qw(gwas_pair_maf_min)),
            gnomad_maf   => pick_field(\%idx, \@f, qw(gnomad_maf)),
            gnomad_pops  => pick_field(\%idx, \@f, qw(gnomad_pops)),
            maf_filter_decision => pick_field(\%idx, \@f, qw(maf_filter_decision)),
            maf_filter_reason => pick_field(\%idx, \@f, qw(maf_filter_reason)),
        };
    }
    close $fh;
    @hits = sort {
        (($a->{hit_order} || 9e9) <=> ($b->{hit_order} || 9e9))
            ||
        chr_order($a->{CHR}) <=> chr_order($b->{CHR})
            ||
        ($a->{BP} || 0) <=> ($b->{BP} || 0)
    } @hits;
    return @hits;
}

sub read_manifest_tsv {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    my %metrics;
    my $header = <$fh>;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my ($metric, $value) = split /\t/, $line, 2;
        $metrics{$metric} = $value;
        $metrics{uc($metric)} = $value;
        $metrics{lc($metric)} = $value;
    }
    close $fh;
    return \%metrics;
}

sub locate_snp_in_gzip {
    my ($path, $target_snp) = @_;
    return unless $path && -s $path && $target_snp;
    my $fh = IO::Uncompress::Gunzip->new($path)
        or die "Cannot read $path: $GunzipError\n";
    my $header = <$fh>;
    return unless defined $header;
    chomp $header;
    my @cols = split /\t/, $header, -1;
    my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
    return unless exists $idx{SNP} && exists $idx{CHR} && exists $idx{BP};
    while (my $line = <$fh>) {
        chomp $line;
        my @f = split /\t/, $line, -1;
        next unless ($f[ $idx{SNP} ] // '') eq $target_snp;
        close $fh;
        return ($f[ $idx{CHR} ], $f[ $idx{BP} ]);
    }
    close $fh;
    return;
}

sub pick_field {
    my ($idx, $fields, @names) = @_;
    for my $name (@names) {
        next unless exists $idx->{$name};
        return $fields->[ $idx->{$name} ];
    }
    return '';
}

sub enrich_hits_from_local_top_hits_csv {
    my (%args) = @_;
    my $csv_name = $args{runner}{LOCAL_TOP_HITS_CSV_BASENAME} || '';
    return unless $csv_name;
    my @candidate_paths = (
        File::Spec->catfile($args{output_dir}, $csv_name),
        File::Spec->catfile($Bin, $csv_name),
        $csv_name,
    );
    my ($path) = grep { defined $_ && -f $_ } @candidate_paths;
    return unless $path;
    open my $fh, '<', $path or return;
    my $header = <$fh>;
    return unless defined $header;
    chomp $header;
    my @cols = split /,/, $header, -1;
    my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
    return unless exists $idx{SNP};
    my %by_snp;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = split /,/, $line, -1;
        my $snp = $f[ $idx{SNP} ] // '';
        next unless length $snp;
        $by_snp{$snp} = {
            chr       => (exists $idx{CHR} ? ($f[ $idx{CHR} ] // '') : ''),
            bp        => (exists $idx{BP} ? ($f[ $idx{BP} ] // '') : ''),
            gene     => (exists $idx{gene} ? ($f[ $idx{gene} ] // '') : ''),
            snp_gene => (exists $idx{snp_gene} ? ($f[ $idx{snp_gene} ] // '') : ''),
            gene_source => (exists $idx{gene_source} ? ($f[ $idx{gene_source} ] // '') : ''),
        };
    }
    close $fh;
    for my $hit (@{ $args{hits} || [] }) {
        next unless $hit->{SNP} && exists $by_snp{ $hit->{SNP} };
        $hit->{CHR} ||= $by_snp{ $hit->{SNP} }{chr};
        $hit->{BP} ||= $by_snp{ $hit->{SNP} }{bp};
        $hit->{gene} ||= $by_snp{ $hit->{SNP} }{gene};
        $hit->{snp_gene} ||= $by_snp{ $hit->{SNP} }{snp_gene};
        $hit->{gene_source} ||= $by_snp{ $hit->{SNP} }{gene_source};
        $hit->{gene} = '' if defined $hit->{gene} && $hit->{gene} =~ /^(?:NA|N\/A|null|\.)$/i;
    }
}

sub ensure_tabix_ready_long_source {
    my (%args) = @_;
    my $source_long = $args{source_long} || '';
    return '' unless $source_long && -s $source_long;
    return $source_long if has_tabix_index_wrapper($source_long);

    my $source_basename = basename($source_long);
    $source_basename =~ s/\.(?:gz|bgz|bgzip)$//i;
    my @indexed_candidates;
    if ($source_basename =~ /\.tsv$/i) {
        (my $legacy_basename = $source_basename) =~ s/\.tsv$//i;
        push @indexed_candidates, File::Spec->catfile(
            $args{output_dir},
            $legacy_basename . '.tabix_ready.tsv.gz',
        );
    }
    my $indexed_local = File::Spec->catfile(
        $args{output_dir},
        $source_basename . '.tabix_ready.tsv.gz',
    );
    push @indexed_candidates, $indexed_local;
    if (!$args{force}) {
        for my $candidate (@indexed_candidates) {
            next unless -s $candidate && has_tabix_index_wrapper($candidate);
            print "[skip] reusing indexed long GWAS source $candidate\n";
            return $candidate;
        }
    }

    my @cmd = (
        $^X,
        File::Spec->catfile($Bin, 'bgzip_tabix_diff_gwas.pl'),
        '--input', windows_to_mnt_path($source_long),
        '--output', windows_to_mnt_path($indexed_local),
        '--seq', 1,
        '--start', 2,
        '--end', 2,
    );
    run_cmd(\@cmd, 'build bgzip/tabix indexed long GWAS source for fast local extraction');
    die "Indexed long GWAS source was not created: $indexed_local\n"
        unless -s $indexed_local && has_tabix_index_wrapper($indexed_local);
    return $indexed_local;
}

sub has_tabix_index_wrapper {
    my ($path) = @_;
    return 0 unless $path && -s $path;
    my $data_mtime = (stat($path))[9] || 0;
    for my $index_path ("$path.tbi", "$path.csi") {
        next unless -s $index_path;
        return 1 if ((stat($index_path))[9] || 0) >= $data_mtime;
    }
    if ($path =~ /\.(?:gz|bgz|bgzip)$/i) {
        (my $stem = $path) =~ s/\.(?:gz|bgz|bgzip)$//i;
        for my $index_path ("$stem.tbi", "$stem.csi") {
            next unless -s $index_path;
            return 1 if ((stat($index_path))[9] || 0) >= $data_mtime;
        }
    }
    return 0;
}

sub extract_gene_from_snp_gene {
    my ($text) = @_;
    return '' unless defined $text && length $text;
    my ($snp, $gene) = split /:/, $text, 2;
    return '' unless defined $gene;
    return '' if $gene =~ /^(?:NA|N\/A|null|\.)$/i;
    return $gene;
}

sub parse_target_snp_gene_overrides {
    my (%args) = @_;
    my %overrides;
    my @sources = grep { defined $_ } ($args{spec_value}, $args{cli_value});
    for my $source (@sources) {
        if (ref($source) eq 'HASH') {
            for my $snp (keys %{$source}) {
                my $gene = trim($source->{$snp});
                next unless length $gene;
                $overrides{ uc(trim($snp)) } = $gene if length trim($snp);
            }
            next;
        }
        next if ref $source;
        for my $pair (split /[,\r\n]+/, $source) {
            $pair = trim($pair);
            next unless length $pair;
            my ($snp, $gene) = split /:/, $pair, 2;
            die "Invalid --target-snp-genes entry '$pair'. Expected SNP:GENE.\n"
                unless defined $snp && defined $gene;
            $snp = trim($snp);
            $gene = trim($gene);
            die "Invalid --target-snp-genes entry '$pair'. Expected SNP:GENE.\n"
                unless length($snp) && length($gene);
            $overrides{ uc($snp) } = $gene;
        }
    }
    return \%overrides;
}

sub apply_target_snp_gene_overrides {
    my (%args) = @_;
    my $overrides = $args{overrides} || {};
    return unless %{$overrides};
    for my $hit (@{ $args{hits} || [] }) {
        my $snp = uc($hit->{SNP} || '');
        next unless length $snp && exists $overrides->{$snp};
        my $gene = $overrides->{$snp};
        next unless defined $gene && length $gene;
        $hit->{gene} = $gene;
        $hit->{snp_gene} = ($hit->{SNP} || '') . ":$gene";
        $hit->{gene_source} = 'USER';
    }
}

sub normalize_requested_plots {
    my ($plots, $step_args) = @_;
    my %requested = (
        plot_manhattan       => 0,
        plot_local_manhattan => 0,
        plot_local_gtf       => 0,
        plot_forest          => 0,
    );
    for my $item (split /,/, ($plots || '')) {
        my $v = trim($item);
        $requested{plot_manhattan} = 1 if $v eq 'manhattan';
        $requested{plot_local_manhattan} = 1 if $v eq 'local_manhattan';
        $requested{plot_local_gtf} = 1 if $v eq 'local_gtf';
        $requested{plot_forest} = 1 if $v eq 'forest';
    }
    for my $step (@{ $step_args || [] }) {
        my $v = trim($step);
        $requested{$v} = 1 if exists $requested{$v};
    }
    return %requested;
}

sub run_cmd {
    my ($cmd, $label) = @_;
    print "[run] $label\n";
    print "  " . join(' ', map { quote_for_log($_) } @{$cmd}) . "\n";
    system @{$cmd};
    my $exit = $? >> 8;
    die "Command failed for $label (exit $exit)\n" if $exit != 0;
}

sub format_elapsed_seconds {
    my ($elapsed) = @_;
    $elapsed = 0 unless defined $elapsed;
    $elapsed = 0 if $elapsed < 0;
    my $hours = int($elapsed / 3600);
    my $minutes = int(($elapsed % 3600) / 60);
    my $seconds = $elapsed - ($hours * 3600) - ($minutes * 60);
    return sprintf('%dh %02dm %.1fs', $hours, $minutes, $seconds) if $hours > 0;
    return sprintf('%dm %.1fs', $minutes, $seconds) if $minutes > 0;
    return sprintf('%.1fs', $seconds);
}

sub quote_for_log {
    my ($text) = @_;
    return '""' unless defined $text && length $text;
    return $text if $text !~ /[\s"]/;
    $text =~ s/"/\\"/g;
    return qq{"$text"};
}

sub load_json {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    local $/;
    my $json = <$fh>;
    close $fh;
    return decode_json($json);
}

sub runner_wide_cache_is_reusable {
    my ($runner) = @_;
    return (0, 'runner config is missing') unless ref($runner) eq 'HASH';

    my $data_path = localize_path($runner->{DATA_GZ} || '');
    return (0, 'DATA_GZ is missing from runner config') unless length $data_path;
    return (0, "wide subset missing: $data_path") unless -s $data_path;

    my $extractor_cfg_path = localize_path($runner->{EXTRACTOR_CONFIG_JSON} || '');
    my $manifest_path = '';
    if ($extractor_cfg_path && -s $extractor_cfg_path) {
        my $extractor_cfg = load_json($extractor_cfg_path);
        $manifest_path = localize_path($extractor_cfg->{manifest} || '');
    }
    $manifest_path ||= derive_manifest_from_wide_output($data_path);
    return (0, "wide subset manifest is missing: $manifest_path") unless length $manifest_path && -s $manifest_path;

    return wide_output_matches_manifest(
        data_path     => $data_path,
        manifest_path => $manifest_path,
    );
}

sub derive_manifest_from_wide_output {
    my ($data_path) = @_;
    return '' unless defined $data_path && length $data_path;
    return $1 . '.manifest.tsv' if $data_path =~ /^(.*)\.tsv\.gz$/i;
    return $1 . '.manifest.tsv' if $data_path =~ /^(.*)\.tsv$/i;
    return '';
}

sub wide_output_matches_manifest {
    my (%args) = @_;
    my $data_path = $args{data_path} || '';
    my $manifest_path = $args{manifest_path} || '';
    return (0, 'wide subset path is missing') unless length $data_path;
    return (0, 'wide manifest path is missing') unless length $manifest_path;

    my %metric;
    open my $mf, '<', $manifest_path or return (0, "cannot read manifest $manifest_path: $!");
    while (my $line = <$mf>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my ($k, $v) = split /\t/, $line, 2;
        next unless defined $k && defined $v;
        $metric{$k} = $v;
    }
    close $mf;

    if (exists $metric{target_snp} || exists $metric{rows_in_window}) {
        my $target = $metric{target_snp} // 'unknown_target';
        return (0, "manifest records a single-SNP/local-window artifact for $target");
    }

    my $expected = $metric{rows_written};
    return (0, "manifest $manifest_path does not record rows_written")
      unless defined $expected && $expected =~ /^\d+$/;

    my $actual = count_data_rows($data_path);
    return (0, "wide subset row count mismatch (actual=$actual, manifest=$expected)")
      if $actual != $expected;

    return (1, '');
}

sub count_data_rows {
    my ($path) = @_;
    my $fh;
    if ($path =~ /\.gz$/i) {
        $fh = IO::Uncompress::Gunzip->new($path)
          or die "Cannot read $path: $GunzipError\n";
    }
    else {
        open $fh, '<', $path or die "Cannot read $path: $!\n";
    }
    my $rows = -1; # subtract the header
    while (my $line = <$fh>) {
        $rows++;
    }
    close $fh;
    $rows = 0 if $rows < 0;
    return $rows;
}

sub find_gnuplot_exe {
    my @candidates = (
        $ENV{GNUPLOT_EXE},
        'gnuplot',
        'C:/Users/cheng/Downloads/strawberry-perl-5.42.0.1-64bit-PDL/c/bin/gnuplot.exe',
        'C:/Users/cheng/Downloads/strawberry-perl-5.32.1.1-64bit-PDL/c/bin/gnuplot.exe',
        'C:/Users/cheng/Downloads/cygwin-portable-20210411/cygwin-portable/App/cygwin/bin/gnuplot-base.exe',
        'C:/Users/cheng/Downloads/cygwin-portable-20210411/cygwin-portable/App/cygwin/bin/gnuplot-X11.exe',
    );
    for my $cand (@candidates) {
        next unless defined $cand && length $cand;
        if ($cand eq 'gnuplot') {
            return $cand if command_exists($cand);
            next;
        }
        return $cand if -f $cand;
    }
    return 'gnuplot';
}

sub command_exists {
    my ($cmd) = @_;
    return 0 unless defined $cmd && length $cmd;
    if ($^O eq 'MSWin32') {
        return system("where $cmd >NUL 2>&1") == 0 ? 1 : 0;
    }
    return system('sh', '-lc', "command -v '$cmd' >/dev/null 2>&1") == 0 ? 1 : 0;
}

sub gunplotize_name {
    my ($name) = @_;
    $name ||= 'gunplot_output';
    $name =~ s/_SAS_/_GUNPLOT_/g;
    if ($name !~ /_GUNPLOT_/ && $name !~ /GUNPLOT/i) {
        my ($stem, $ext) = $name =~ /^(.*?)(\.[^.]+)?$/;
        $name = $stem . '_GUNPLOT' . ($ext || '');
    }
    return $name;
}

sub gunplot_title {
    my ($title) = @_;
    $title ||= 'Plot output';
    return $title;
}

sub localize_path {
    my ($path) = @_;
    return '' unless defined $path && length $path;
    $path =~ s{\\}{/}g;
    if ($^O =~ /^(?:cygwin|MSWin32)$/i) {
        if ($path =~ m{^/mnt/([a-zA-Z])/(.*)$}) {
            return uc($1) . ':/' . $2;
        }
        if ($path =~ m{^/cygdrive/([a-zA-Z])/(.*)$}) {
            return uc($1) . ':/' . $2;
        }
        return $path;
    }
    if ($path =~ m{^([a-zA-Z]):/(.*)$}) {
        my ($drive, $rest) = (lc($1), $2);
        $rest =~ s{^/+}{};
        return "/mnt/$drive/$rest";
    }
    return $path;
}

sub resolve_portable_workdir {
    my ($configured, $fallback) = @_;
    $fallback = localize_path($fallback || $Bin);
    return localize_path($ENV{PIPELINE_WORKDIR})
        if defined($ENV{PIPELINE_WORKDIR}) && length($ENV{PIPELINE_WORKDIR});
    my $resolved = localize_path($configured || $fallback);
    return $resolved if $ENV{PIPELINE_RESPECT_SPEC_WORKDIR};
    if (length($resolved) && length($fallback) && $resolved ne $fallback) {
        warn "[info] Overriding spec workdir '$resolved' with current pipeline directory '$fallback' for portable execution.\n";
        return $fallback;
    }
    return $resolved;
}

sub windows_to_mnt_path {
    my ($path) = @_;
    return '' unless defined $path && length $path;
    $path =~ s{\\}{/}g;
    return $path if $path =~ m{^/mnt/[A-Za-z]/};
    return $path if $path =~ m{^/cygdrive/[A-Za-z]/};
    if ($path =~ m{^([A-Za-z]):/(.*)$}) {
        my ($drive, $rest) = (lc($1), $2);
        $rest =~ s{^/+}{};
        return "/mnt/$drive/$rest";
    }
    return $path;
}

sub unix_join {
    return join('/', map { my $v = $_; $v =~ s{/$}{}; $v } @_);
}

sub safe_name {
    my ($text) = @_;
    $text //= '';
    $text =~ s/[^A-Za-z0-9._-]+/_/g;
    $text =~ s/^_+|_+$//g;
    return length($text) ? $text : 'item';
}

sub compute_gtf_height {
    my ($track_count) = @_;
    $track_count ||= 1;
    return 900 + ($track_count * 55);
}

sub repeated_panel_ytics_wrapper {
    my ($track_count, $top_logp) = @_;
    my @ticks;
    for my $i (0 .. ($track_count - 1)) {
        my $base = $i * $top_logp;
        for my $tick (0, 2, 4, 6, 8) {
            next if $tick > $top_logp;
            push @ticks, sprintf('"%d" %.3f', $tick, $base + $tick);
        }
    }
    return @ticks;
}

sub safe_neglog10_text {
    my ($p_text) = @_;
    my $p = $p_text;
    return 0 unless defined $p && $p =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/ && $p > 0;
    return -log($p) / log(10);
}

sub chromosome_color_wrapper {
    my ($chr) = @_;
    my @palette = sas_chr_palette_wrapper();
    my $ord = chr_palette_index_wrapper($chr);
    return $palette[($ord - 1) % @palette];
}

sub sas_chr_palette_wrapper {
    return (
        '#0072bd', '#d95319', '#edb120', '#7e2f8e',
        '#77ac30', '#4dbeee', '#a2142f',
    );
}

sub chr_palette_index_wrapper {
    my ($chr) = @_;
    return 23 if defined $chr && $chr =~ /^(?:X|23)$/i;
    return $chr if defined $chr && $chr =~ /^\d+$/ && $chr > 0;
    return 1;
}

sub nearest_gene_label_for_item {
    my ($item) = @_;
    my $gene = $item->{gene} || '';
    $gene = extract_gene_from_snp_gene($item->{snp_gene}) if !length($gene);
    return $gene if length $gene;

    my $manifest = read_manifest_tsv($item->{manifest});
    my $bp = $item->{bp} || $manifest->{TARGET_BP} || $manifest->{target_bp} || $manifest->{BP} || $manifest->{bp};
    my $chr = $item->{chr} || $manifest->{TARGET_CHR} || $manifest->{target_chr} || $manifest->{CHR} || $manifest->{chr};
    return '' unless defined $bp;

    my $genes_tsv = '';
    if ($item->{image}) {
        $genes_tsv = $item->{image};
        $genes_tsv =~ s/\.png$/.genes.tsv/i;
    }
    if ($genes_tsv && -f $genes_tsv) {
        open my $gfh, '<', $genes_tsv or return '';
        my $gheader = <$gfh>;
        if (defined $gheader) {
            my ($best_gene, $best_dist);
            while (my $line = <$gfh>) {
                chomp $line;
                $line =~ s/\r$//;
                next unless length $line;
                my ($type, $gene_name, $start, $end) = split /\t/, $line, -1;
                next unless defined $type && lc($type) eq 'gene';
                next unless defined $gene_name && length $gene_name;
                next unless defined $start && defined $end && $start =~ /^\d+$/ && $end =~ /^\d+$/;
                my $dist = ($bp >= $start && $bp <= $end) ? 0 : ($bp < $start ? $start - $bp : $bp - $end);
                if (!defined $best_dist || $dist < $best_dist || ($dist == $best_dist && $gene_name cmp ($best_gene || '') < 0)) {
                    $best_dist = $dist;
                    $best_gene = $gene_name;
                }
            }
            close $gfh;
            return $best_gene if defined $best_gene && length $best_gene;
        } else {
            close $gfh;
        }
    }

    return '' unless defined $chr;

    my $gtf_file = $item->{gtf_file};
    if (!$gtf_file || !-f $gtf_file) {
        my $guess = $item->{image};
        $guess =~ s/\.png$/.gtf.tsv/i if defined $guess;
        $gtf_file = $guess if defined $guess && -f $guess;
    }
    return nearest_gene_from_gtf_file(
        gtf_file => $gtf_file,
        chr      => $chr,
        bp       => $bp,
    );
}

sub nearest_gene_from_gtf_file {
    my (%args) = @_;
    my $gtf_file = $args{gtf_file} || '';
    my $chr = $args{chr};
    my $bp = $args{bp};
    return '' unless $gtf_file && -f $gtf_file && defined $bp;
    my $query_chr = normalize_chr_label($chr);
    my $restrict_chr = length $query_chr ? 1 : 0;

    open my $fh, '<', $gtf_file or return '';
    my $header = <$fh>;
    return '' unless defined $header;
    chomp $header;
    my @cols = split /\t/, $header, -1;
    my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
    for my $need (qw(chr genesymbol gene start end type)) {
        return '' unless exists $idx{$need};
    }
    my ($best_gene, $best_dist);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        next unless lc(($f[ $idx{type} ] // '')) eq 'gene';
        next if $restrict_chr && normalize_chr_label($f[ $idx{chr} ] // '') ne $query_chr;
        my $gene_name = $f[ $idx{genesymbol} ] || $f[ $idx{gene} ] || '';
        next unless length $gene_name;
        my $start = $f[ $idx{start} ];
        my $end = $f[ $idx{end} ];
        next unless defined $start && defined $end && $start =~ /^\d+$/ && $end =~ /^\d+$/;
        my $dist = ($bp >= $start && $bp <= $end) ? 0 : ($bp < $start ? $start - $bp : $bp - $end);
        if (!defined $best_dist || $dist < $best_dist || ($dist == $best_dist && $gene_name cmp $best_gene < 0)) {
            $best_dist = $dist;
            $best_gene = $gene_name;
        }
    }
    close $fh;
    return $best_gene || '';
}

sub normalize_chr_label {
    my ($chr) = @_;
    $chr = '' unless defined $chr;
    $chr = trim($chr);
    $chr =~ s/^chr//i;
    return uc($chr);
}

sub infer_cached_locus_gtf_path {
    my (%args) = @_;
    my $snp = $args{snp} || '';
    return '' unless length $snp;
    my $safe_snp = safe_name($snp);
    my $safe_window = safe_name($args{window_bp});
    my $npc_flag = ($args{runner} && $args{runner}{LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}) ? 1 : 0;
    my $gtf_path = File::Spec->catfile(
        $args{output_dir},
        "gunplot_locus_${safe_snp}_window_${safe_window}_npc${npc_flag}.gtf.tsv",
    );
    return $gtf_path if -f $gtf_path;
    my ($fallback) = glob(File::Spec->catfile($args{output_dir}, "gunplot_locus_${safe_snp}_window_*_npc${npc_flag}.gtf.tsv"));
    return ($fallback && -f $fallback) ? $fallback : '';
}

sub nearest_gene_for_snp_from_cached_gtf {
    my (%args) = @_;
    my $gtf_file = infer_cached_locus_gtf_path(
        output_dir => $args{output_dir},
        snp        => $args{snp},
        window_bp  => $args{window_bp},
        runner     => {},
    );
    return '' unless $gtf_file;
    return nearest_gene_from_gtf_file(
        gtf_file => $gtf_file,
        chr      => '',
        bp       => $args{bp},
    );
}

sub ensure_genes_tsv_from_gtf {
    my (%args) = @_;
    my $gtf_file = $args{gtf_file} || '';
    my $gene_tsv = $args{gene_tsv} || '';
    return unless $gtf_file && -f $gtf_file && $gene_tsv;
    return if -s $gene_tsv;

    open my $in, '<', $gtf_file or return;
    my $header = <$in>;
    return unless defined $header;
    chomp $header;
    my @cols = split /\t/, $header, -1;
    my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
    for my $need (qw(type start end genesymbol gene)) {
        return unless exists $idx{$need};
    }

    open my $out, '>', $gene_tsv or return;
    print {$out} join("\t", qw(TYPE GENE START END LANE)), "\n";
    my %seen;
    while (my $line = <$in>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        next unless lc(($f[$idx{type}] // '')) eq 'gene';
        my $gene = $f[$idx{genesymbol}] || $f[$idx{gene}] || '';
        next unless length $gene;
        my $start = $f[$idx{start}] // '';
        my $end = $f[$idx{end}] // '';
        next unless $start =~ /^\d+$/ && $end =~ /^\d+$/;
        my $key = join("\t", $gene, $start, $end);
        next if $seen{$key}++;
        print {$out} join("\t", 'gene', $gene, $start, $end, 0), "\n";
    }
    close $out;
    close $in;
}

sub mini_gtf_genes_for_column {
    my ($item) = @_;
    return ([], 0) unless $item && $item->{gtf_file} && -f $item->{gtf_file};
    my $manifest = read_manifest_tsv($item->{manifest});
    my $bp = $item->{bp} || $manifest->{TARGET_BP} || $manifest->{target_bp} || $manifest->{BP} || $manifest->{bp};
    my $chr = $item->{chr} || $manifest->{TARGET_CHR} || $manifest->{target_chr} || $manifest->{CHR} || $manifest->{chr};
    my $window_bp = $manifest->{WINDOW_BP} || $manifest->{window_bp} || 0;
    my $start = $manifest->{WINDOW_START} || $manifest->{window_start} || ($bp - $window_bp);
    my $end   = $manifest->{WINDOW_END} || $manifest->{window_end} || ($bp + $window_bp);
    return ([], 0) unless defined $start && defined $end && $end > $start;

    open my $fh, '<', $item->{gtf_file} or return ([], 0);
    my $header = <$fh>;
    return ([], 0) unless defined $header;
    chomp $header;
    my @cols = split /\t/, $header, -1;
    my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
    for my $need (qw(chr genesymbol gene start end type)) {
        return ([], 0) unless exists $idx{$need};
    }

    my @genes;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        next unless lc(($f[ $idx{type} ] // '')) eq 'gene';
        next unless ($f[ $idx{chr} ] // '') eq $chr;
        my $gene_name = $f[ $idx{genesymbol} ] || $f[ $idx{gene} ] || '';
        next unless length $gene_name;
        my $gs = $f[ $idx{start} ];
        my $ge = $f[ $idx{end} ];
        next unless defined $gs && defined $ge && $gs =~ /^\d+$/ && $ge =~ /^\d+$/;
        next if $ge < $start || $gs > $end;
        my $mid = ($gs + $ge) / 2;
        my $dist = defined $bp ? abs($mid - $bp) : 0;
        push @genes, {
            gene => $gene_name,
            start => $gs,
            end => $ge,
            dist => $dist,
        };
    }
    close $fh;

    @genes = sort {
        $a->{dist} <=> $b->{dist}
            ||
        (($b->{end} - $b->{start}) <=> ($a->{end} - $a->{start}))
            ||
        ($a->{gene} cmp $b->{gene})
    } @genes;

    my %seen;
    @genes = grep { !$seen{ $_->{gene} }++ } @genes;
    splice(@genes, 40) if @genes > 40;
    @genes = sort { $a->{start} <=> $b->{start} || $a->{end} <=> $b->{end} } @genes;

    my @palette = (
        '#f4a6a6', '#f7c97f', '#d8d35f', '#9fd27a', '#71c9b8',
        '#86b6f6', '#a995e8', '#de9de6', '#f4a8c4', '#c8b08b',
        '#ef8b62', '#9ec3a5',
    );
    my @lane_ends;
    my $lane_count = 0;
    for my $i (0 .. $#genes) {
        my $g = $genes[$i];
        my $slot = $item->{slot_index} // 0;
        my ($col_min, $col_max) = locus_column_bounds($slot);
        my $x1 = locus_bp_to_column_x($slot, $start, $end, $g->{start});
        my $x2 = locus_bp_to_column_x($slot, $start, $end, $g->{end});
        my $min_w = 0.028;
        if (($x2 - $x1) < $min_w) {
            my $mid = ($x1 + $x2) / 2;
            $x1 = $mid - $min_w / 2;
            $x2 = $mid + $min_w / 2;
        }
        $x1 = $col_min if $x1 < $col_min;
        $x2 = $col_max if $x2 > $col_max;
        if (($x2 - $x1) < $min_w) {
            my $mid = ($x1 + $x2) / 2;
            $x1 = $mid - $min_w / 2;
            $x2 = $mid + $min_w / 2;
            $x1 = $col_min if $x1 < $col_min;
            $x2 = $col_max if $x2 > $col_max;
        }
        my $lane = allocate_compact_lane(\@lane_ends, $x1, $x2, 0.018);
        $lane_count = $lane + 1 if ($lane + 1) > $lane_count;
        $g->{lane} = $lane;
        $g->{x1} = $x1;
        $g->{x2} = $x2;
        $g->{xc} = ($x1 + $x2) / 2;
        $g->{color} = $palette[$i % @palette];
    }
    return (\@genes, $lane_count);
}

sub allocate_compact_lane {
    my ($lane_ends, $start, $end, $pad) = @_;
    $pad ||= 0;
    for my $i (0 .. $#{$lane_ends}) {
        next if $start <= ($lane_ends->[$i] + $pad);
        $lane_ends->[$i] = $end;
        return $i;
    }
    push @{$lane_ends}, $end;
    return $#{$lane_ends};
}

sub locus_column_bounds {
    my ($slot_index) = @_;
    $slot_index ||= 0;
    return ($slot_index + 0.08, $slot_index + 0.92);
}

sub locus_bp_to_column_x {
    my ($slot_index, $window_start, $window_end, $bp) = @_;
    my ($col_min, $col_max) = locus_column_bounds($slot_index);
    return ($col_min + $col_max) / 2
        unless defined $window_start && defined $window_end && $window_end > $window_start;
    my $frac = ($bp - $window_start) / ($window_end - $window_start);
    $frac = 0 if $frac < 0;
    $frac = 1 if $frac > 1;
    return $col_min + ($col_max - $col_min) * $frac;
}

sub normalize_local_manhattan_annotation_mode {
    my ($mode, $n_cols) = @_;
    $mode = lc(trim($mode || 'gtf'));
    $mode = 'labels' unless $mode =~ /^(?:labels|gtf|auto|none)$/;
    return 'none' if $mode eq 'none';
    if ($mode eq 'auto') {
        return ($n_cols && $n_cols < 5) ? 'gtf' : 'labels';
    }
    return $mode;
}

sub html_escape {
    my ($text) = @_;
    $text //= '';
    $text =~ s/&/&amp;/g;
    $text =~ s/</&lt;/g;
    $text =~ s/>/&gt;/g;
    $text =~ s/"/&quot;/g;
    return $text;
}

sub resolve_runner_differential_threshold_ladder {
    my ($runner) = @_;
    $runner ||= {};
    my $primary = trim($runner->{TOP_HIT_SIGNAL_THRSHD} || '1e-6');
    my $runner_ladder = trim($runner->{TOP_HIT_SIGNAL_THRSHDS} || '');
    my $fallback = trim($runner->{TOP_HIT_SIGNAL_THRSHD_FALLBACK} || '');
    my @vals = grep { length } split /[,\s]+/, $runner_ladder;
    @vals = ($primary) unless @vals;
    if (@vals <= 1) {
        $fallback = '1e-5'
            if !length($fallback)
            && normalized_numeric_text($primary) eq normalized_numeric_text('1e-6');
        push @vals, $fallback if length($fallback) && !grep { $_ eq $fallback } @vals;
    }
    return join(' ', unique_threshold_values(@vals));
}

sub unique_threshold_values {
    my @values = @_;
    my @uniq;
    for my $value (@values) {
        next unless defined $value;
        $value = trim($value);
        next unless length $value;
        push @uniq, $value unless grep { $_ eq $value } @uniq;
    }
    return @uniq;
}

sub normalized_numeric_text {
    my ($value) = @_;
    return '' unless defined $value;
    $value = trim($value);
    return '' unless length $value;
    my $num = eval { 0 + $value };
    return $value if $@;
    return sprintf('%.12g', $num);
}

sub escape_gp {
    my ($text) = @_;
    $text //= '';
    $text =~ s/\\/\\\\/g;
    $text =~ s/'/\\'/g;
    $text =~ s/"/\\"/g;
    return $text;
}

sub italic_font_spec_gp {
    my ($size) = @_;
    $size ||= 11;
    my @candidates = (
        'C:/Windows/Fonts/timesi.ttf',
        'C:/Windows/Fonts/ariali.ttf',
    );
    for my $cand (@candidates) {
        return $cand . "," . $size if -f $cand;
    }
    return "," . $size;
}

sub csv_escape {
    my ($text) = @_;
    $text //= '';
    if ($text =~ /[",\n]/) {
        $text =~ s/"/""/g;
        return qq{"$text"};
    }
    return $text;
}

sub trim {
    my ($x) = @_;
    $x //= '';
    $x =~ s/^\s+//;
    $x =~ s/\s+$//;
    return $x;
}

sub chr_order {
    my ($chr) = @_;
    return 23 if defined $chr && uc($chr) eq 'X';
    return $chr if defined $chr && $chr =~ /^\d+$/;
    return 10_000;
}

sub is_chr_x {
    my ($chr) = @_;
    return 0 unless defined $chr;
    $chr =~ s/^\s+|\s+$//g;
    $chr =~ s/^chr//i;
    return $chr =~ /^(?:X|23)$/i ? 1 : 0;
}
