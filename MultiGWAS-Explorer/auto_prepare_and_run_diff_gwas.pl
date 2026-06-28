#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
BEGIN {
    my $self_dir = $FindBin::RealBin || $FindBin::Bin || '.';
    $self_dir =~ s{\\}{/}g;
    require lib;
    lib->import("$self_dir/DiffGWASDeps");
}
use Getopt::Long qw(GetOptions);
use JSON::PP qw(decode_json encode_json);
use IO::Uncompress::Gunzip qw($GunzipError);
use Digest::MD5 qw(md5_hex);
use Cwd qw(abs_path);
use File::Spec;
use File::Basename qw(basename);
use File::Temp qw(tempfile);
use POSIX qw(strftime);
use Time::HiRes qw(time);
use DiffGWASRawSchema qw(
  normalize_header_name
  resolve_raw_header_aliases
);
use GenomeBuildProfile qw(
  canonicalize_reference_build
  detect_reference_build_profile
);

BEGIN {
    my $old = select STDERR;
    $| = 1;
    select $old;
}

my $spec_file = '';
my $gwas_dir = '';
my $spec_out = '';
my $raw_column_alias_config = '';
my $mode = 'full';
my $skip_plots = 0;
my $plots = 'manhattan,local_manhattan,local_gtf';
my $force = 0;
my @step_args;
my $from_step = '';
my $to_step = '';
my $list_steps = 0;
my $print_spec_example = 0;
my $print_columns_help = 0;
my $generate_spec_only = 0;
my $preview_spec = 0;
my $project_tag_override = '';
my $artifact_stem_override = '';
my $reference_build_override = '';
my $top_hit_dist_bp_override = '';
my $top_hit_max_loci_override;
my $local_max_hits_per_fig_override = 0;
my $local_gtf_window_bp_override = '';
my $local_manhattan_angle4xaxis_label_override = '';
my $local_manhattan_xgrp_y_pos_override = '';
my $local_manhattan_yoffset_top_override = '';
my $local_manhattan_yoffset_bottom_override = '';
my $local_manhattan_fontsize_override = '';
my $local_manhattan_y_axis_label_size_override = '';
my $local_manhattan_y_axis_value_size_override = '';
my $target_snps_override = '';
my $target_snp_genes_override = '';
my $display_gwas_override = '';
my $sas_oda_account_override = '';
my $sas_oda_password_override = '';
my $prompt_sas_oda_auth_override = 0;
my $emit_local_sas_scripts = 0;
my $local_sas_only = 0;
my $local_gtf_label_snps_override = '';
my $local_gtf_label_layout_override = '';
my $local_gtf_yaxis_offset4max_override = '';
my $local_gtf_yoffset4textlabels_override = '';
my $exclude_non_protein_coding_genes_in_local_gtf = 0;
my $get_common_associations;
my $common_assoc_top_hit_threshold_override = '';
my $EXTERNAL_GZIP = '';
my %step_flag = (
    merge_raw                 => 0,
    sort_long                 => 0,
    diff_pairs                => 0,
    standardize_diff          => 0,
    extract_wide_subset       => 0,
    plot_manhattan            => 0,
    plot_local_manhattan      => 0,
    plot_local_gtf            => 0,
    plot_forest               => 0,
    cleanup_shared_plot_data  => 0,
);

if (!@ARGV) {
    print full_help();
    exit 0;
}

for my $arg (@ARGV) {
    next unless defined $arg;
    next if $arg eq '--target-snp-genes';
    next if $arg !~ /^--/;
    if ($arg =~ /target-snp-genes/i || $arg =~ /snp-genes/i) {
        die "Error: Unrecognized option '$arg'. Did you mean '--target-snp-genes'?\n";
    }
}

GetOptions(
    'spec=s'              => \$spec_file,
    'gwas-dir=s'          => \$gwas_dir,
    'spec-out=s'          => \$spec_out,
    'raw-column-alias-config=s' => \$raw_column_alias_config,
    'mode=s'              => \$mode,
    'skip-plots!'         => \$skip_plots,
    'plots=s'             => \$plots,
    'force!'              => \$force,
    'step=s@'             => \@step_args,
    'from-step=s'         => \$from_step,
    'to-step=s'           => \$to_step,
    'list-steps!'         => \$list_steps,
    'print-spec-example!' => \$print_spec_example,
    'print-columns-help!' => \$print_columns_help,
    'generate-spec-only!' => \$generate_spec_only,
    'preview-spec!'       => \$preview_spec,
    'project-tag=s'       => \$project_tag_override,
    'artifact-stem=s'     => \$artifact_stem_override,
    'reference-build=s'   => \$reference_build_override,
    'top-hit-dist-bp=s'   => \$top_hit_dist_bp_override,
    'top-hit-max-loci=i'  => \$top_hit_max_loci_override,
    'local-max-hits-per-fig=i' => \$local_max_hits_per_fig_override,
    'local-gtf-window-bp=s' => \$local_gtf_window_bp_override,
    'local-manhattan-angle4xaxis-label=s' => \$local_manhattan_angle4xaxis_label_override,
    'local-manhattan-xgrp-y-pos=s' => \$local_manhattan_xgrp_y_pos_override,
    'local-manhattan-yoffset-top=s' => \$local_manhattan_yoffset_top_override,
    'local-manhattan-yoffset-bottom=s' => \$local_manhattan_yoffset_bottom_override,
    'local-manhattan-fontsize=s' => \$local_manhattan_fontsize_override,
    'local-manhattan-y-axis-label-size=s' => \$local_manhattan_y_axis_label_size_override,
    'local-manhattan-y-axis-value-size=s' => \$local_manhattan_y_axis_value_size_override,
    'target-snps=s' => \$target_snps_override,
    'target-snp-genes=s' => \$target_snp_genes_override,
    'display-gwas|display-tracks=s' => \$display_gwas_override,
    'sas-oda-account=s' => \$sas_oda_account_override,
    'sas-oda-password=s' => \$sas_oda_password_override,
    'prompt-sas-oda-auth!' => \$prompt_sas_oda_auth_override,
    'emit-local-sas-scripts!' => \$emit_local_sas_scripts,
    'local-sas-only!' => \$local_sas_only,
    'local-gtf-label-snps=s' => \$local_gtf_label_snps_override,
    'local-gtf-label-layout=s' => \$local_gtf_label_layout_override,
    'local-gtf-yaxis-offset4max=s' => \$local_gtf_yaxis_offset4max_override,
    'local-gtf-yoffset4textlabels=s' => \$local_gtf_yoffset4textlabels_override,
    'exclude-non-protein-coding-genes-in-local-gtf!' => \$exclude_non_protein_coding_genes_in_local_gtf,
    'get-common-associations:s' => \$get_common_associations,
    'common-association-top-hit-threshold=s' => \$common_assoc_top_hit_threshold_override,
    'merge-raw!'                => \$step_flag{merge_raw},
    'sort-long!'                => \$step_flag{sort_long},
    'diff-pairs!'               => \$step_flag{diff_pairs},
    'standardize-diff!'         => \$step_flag{standardize_diff},
    'extract-wide-subset!'      => \$step_flag{extract_wide_subset},
    'plot-manhattan!'           => \$step_flag{plot_manhattan},
    'plot-local-manhattan!'     => \$step_flag{plot_local_manhattan},
    'plot-local-gtf!'           => \$step_flag{plot_local_gtf},
    'plot-forest!'              => \$step_flag{plot_forest},
    'cleanup-shared-plot-data!' => \$step_flag{cleanup_shared_plot_data},
) or die usage();

if ($print_spec_example || $print_columns_help) {
    print full_help(
        show_example => $print_spec_example,
        show_columns => $print_columns_help,
    );
    exit 0;
}

$ENV{PIPELINE_SAS_ODA_ACCOUNT} = $sas_oda_account_override
  if defined $sas_oda_account_override && length $sas_oda_account_override;
$ENV{PIPELINE_SAS_ODA_PASSWORD} = $sas_oda_password_override
  if defined $sas_oda_password_override && length $sas_oda_password_override;
$ENV{PIPELINE_FORCE_SAS_ODA_AUTH_PROMPT} = 1 if $prompt_sas_oda_auth_override;

my $cli_raw_column_aliases = load_alias_override_file($raw_column_alias_config);

if (!length $spec_file && length $gwas_dir) {
    my $draft_spec;my $draft_path;
    ($draft_spec, $draft_path) = infer_spec_from_gwas_dir(
        gwas_dir              => $gwas_dir,
        spec_out              => $spec_out,
        workdir               => script_root_dir(),
        project_tag_override  => $project_tag_override,
        artifact_stem_override => $artifact_stem_override,
        raw_column_aliases    => $cli_raw_column_aliases,
    );
    print STDERR "WARNING: The draft_path is undefined but the inferred draft spec is available in memory.\n"
      unless defined($draft_path) && length($draft_path) > 0;
    if ($preview_spec) {
        print pretty_json($draft_spec), "\n";
        exit 0 if $generate_spec_only || !length $spec_out;
    }
    write_json_if_defined($draft_path, $draft_spec);
    if (defined $draft_path){
        print "Generated comparison spec: $draft_path\n";
    }
    print "Detected groups: " . join(', ', map { $_->{tag} } @{ $draft_spec->{groups} || [] }) . "\n";
    print "Detected pairs: " . join(', ', map { $_->{prefix} } @{ $draft_spec->{pairs} || [] }) . "\n" if ref($draft_spec->{pairs}) eq 'ARRAY';
    if ($generate_spec_only) {
        exit 0;
    }
    $spec_file = $draft_path;
}

die "--spec is required (or provide --gwas-dir to generate one)\n" unless length $spec_file;
my $spec = load_json($spec_file);
$spec->{raw_column_aliases} = merge_alias_override_specs($spec->{raw_column_aliases}, $cli_raw_column_aliases)
  if ref($cli_raw_column_aliases) eq 'HASH' && keys %{$cli_raw_column_aliases};
$spec->{reference_build} = $reference_build_override
  if defined $reference_build_override && length $reference_build_override;
$spec->{local_gtf_window_bp} = $local_gtf_window_bp_override
  if defined $local_gtf_window_bp_override && length $local_gtf_window_bp_override;
$spec->{gtf_yaxis_offset4max} = $local_gtf_yaxis_offset4max_override
  if defined $local_gtf_yaxis_offset4max_override && length $local_gtf_yaxis_offset4max_override;
$spec->{gtf_yoffset4textlabels} = $local_gtf_yoffset4textlabels_override
  if defined $local_gtf_yoffset4textlabels_override && length $local_gtf_yoffset4textlabels_override;
$spec->{include_non_protein_coding_genes_in_local_gtf} = 0
  if $exclude_non_protein_coding_genes_in_local_gtf;
my $workdir = resolve_pipeline_workdir(
    spec     => $spec,
    fallback => script_root_dir(),
);
my $configs_dir = normalize_unix_path(cfg_or($spec, 'configs_dir', "$workdir/configs"));
my $runner_session = cfg_or($spec, 'session_id', 'mysession');
my $bash_path = normalize_unix_path(cfg_or(
    $spec,
    'cygwin_bash',
    '/bin/bash'
));

validate_spec($spec);

my $artifact_stem = cfg_or($spec, 'artifact_stem', cfg_or($spec, 'project_tag', 'diff_gwas'));
my $project_tag = cfg_or($spec, 'project_tag', $artifact_stem);
my $input_dir = normalize_unix_path(cfg_or($spec, 'input_dir', ''));
my $output_dir = normalize_unix_path(cfg_or($spec, 'output_dir', $input_dir || $workdir));
my $source_mode = cfg_or($spec, 'source_mode', 'raw_pgc_vcf_sumstats');
my $reference_build_profile = resolve_reference_build_profile_for_spec(
    spec               => $spec,
    reference_override => $reference_build_override,
    source_mode        => $source_mode,
);
my $threshold = cfg_or($spec, 'threshold', 0.05);
my $window_bp = cfg_or($spec, 'window_bp', 10_000_000);
my $top_hit_threshold = cfg_or($spec, 'top_hit_signal_thrshd', '1e-6');
my $top_hit_threshold_fallback = cfg_or($spec, 'top_hit_signal_thrshd_fallback', '');
my $top_hit_dist_bp = length($top_hit_dist_bp_override)
    ? $top_hit_dist_bp_override
    : cfg_or($spec, 'top_hit_dist_bp', '1e6');
my $top_hit_max_loci = defined($top_hit_max_loci_override)
    ? $top_hit_max_loci_override
    : cfg_or($spec, 'top_hit_max_loci', 0);
my $local_window_bp = cfg_or($spec, 'local_window_bp', '1e7');
my $local_gtf_window_bp = cfg_or($spec, 'local_gtf_window_bp', $local_window_bp);
my $open_result = cfg_or($spec, 'open_result', 0) ? 1 : 0;
my $clean_oda_input = cfg_or($spec, 'clean_oda_input', 1) ? 1 : 0;
my $keep_remote_plot_data = cfg_or($spec, 'keep_remote_plot_data', 0) ? 1 : 0;
my $env_open_result = exists $ENV{OPEN_RESULT} ? ($ENV{OPEN_RESULT} ? 1 : 0) : undef;
my $env_clean_oda_input = exists $ENV{CLEAN_ODA_INPUT} ? ($ENV{CLEAN_ODA_INPUT} ? 1 : 0) : undef;
my $env_keep_remote_plot_data = exists $ENV{KEEP_REMOTE_PLOT_DATA} ? ($ENV{KEEP_REMOTE_PLOT_DATA} ? 1 : 0) : undef;
my $env_skip_data_upload = exists $ENV{SKIP_DATA_UPLOAD} ? ($ENV{SKIP_DATA_UPLOAD} ? 1 : 0) : undef;
$open_result = $env_open_result if defined $env_open_result;
$clean_oda_input = $env_clean_oda_input if defined $env_clean_oda_input;
$keep_remote_plot_data = $env_keep_remote_plot_data if defined $env_keep_remote_plot_data;
if ($local_sas_only) {
    $emit_local_sas_scripts = 1;
}

my $generated = build_generated_paths(
    configs_dir   => $configs_dir,
    output_dir    => $output_dir,
    artifact_stem => $artifact_stem,
);
if ($source_mode eq 'merged_gwas_table') {
    $generated->{wide_output} = "$output_dir/" . safe_name($artifact_stem) . ".merged_plotwide.tsv.gz";
    $generated->{wide_manifest} = "$output_dir/" . safe_name($artifact_stem) . ".merged_plotwide.manifest.tsv";
    $generated->{stdized_output} = '';
    $generated->{stdized_manifest} = '';
}
my $deps_dir = normalize_unix_path(File::Spec->catdir($Bin, 'DiffGWASDeps'));
verify_diff_gwas_deps($deps_dir);
my $oda_helper_unix = resolve_oda_helper_unix($Bin);

my $pair_info = build_pair_info($spec->{pairs});
my @prefixes = @{ $pair_info->{prefix_order} };
my @labels = @{ $pair_info->{labels} };
my $focus_prefix = cfg_or($spec, 'top_hit_focus_prefix', $prefixes[0]);
die "top_hit_focus_prefix $focus_prefix is not one of: " . join(', ', @prefixes) . "\n"
  unless grep { $_ eq $focus_prefix } @prefixes;

my $merge_cfg = build_merge_config($spec, $generated) if $source_mode eq 'raw_pgc_vcf_sumstats';
my $diff_cfg = build_diff_config($spec, $generated) if $source_mode eq 'raw_pgc_vcf_sumstats';
my $preset_cfg = build_preset_config($spec, $generated, $pair_info, $threshold, $window_bp);
my $runner_cfg = build_runner_config(
    spec               => $spec,
    generated          => $generated,
    pair_info          => $pair_info,
    project_tag        => $project_tag,
    reference_build_profile => $reference_build_profile,
    focus_prefix       => $focus_prefix,
    top_hit_threshold  => $top_hit_threshold,
    top_hit_dist_bp    => $top_hit_dist_bp,
    top_hit_max_loci   => $top_hit_max_loci,
    local_window_bp    => $local_window_bp,
    local_gtf_window_bp => $local_gtf_window_bp,
    local_max_hits_per_fig_override => $local_max_hits_per_fig_override,
    local_manhattan_angle4xaxis_label_override => $local_manhattan_angle4xaxis_label_override,
    local_manhattan_xgrp_y_pos_override => $local_manhattan_xgrp_y_pos_override,
    local_manhattan_yoffset_top_override => $local_manhattan_yoffset_top_override,
    local_manhattan_yoffset_bottom_override => $local_manhattan_yoffset_bottom_override,
    local_manhattan_fontsize_override => $local_manhattan_fontsize_override,
    local_manhattan_y_axis_label_size_override => $local_manhattan_y_axis_label_size_override,
    local_manhattan_y_axis_value_size_override => $local_manhattan_y_axis_value_size_override,
    target_snps_override => $target_snps_override,
    target_snp_genes_override => $target_snp_genes_override,
    display_gwas_override => $display_gwas_override,
    local_gtf_label_snps_override => $local_gtf_label_snps_override,
    local_gtf_label_layout_override => $local_gtf_label_layout_override,
    local_gtf_yaxis_offset4max_override => $local_gtf_yaxis_offset4max_override,
    local_gtf_yoffset4textlabels_override => $local_gtf_yoffset4textlabels_override,
    get_common_associations => ((defined $get_common_associations) || cfg_or($spec, 'get_common_associations', 0)),
    common_assoc_top_hit_threshold_override => (
        length($common_assoc_top_hit_threshold_override)
        ? $common_assoc_top_hit_threshold_override
        : (defined($get_common_associations) && length($get_common_associations) ? $get_common_associations : '')
    ),
);

write_json_if_defined($generated->{merge_config}, $merge_cfg) if $merge_cfg;
write_json_if_defined($generated->{diff_config}, $diff_cfg) if $diff_cfg;
write_json_if_defined($generated->{preset_config}, $preset_cfg);
write_json_if_defined($generated->{runner_config}, $runner_cfg);

my %summary = (
    spec_file => $spec_file,
    mode => $mode,
    source_mode => $source_mode,
    reference_build => $reference_build_profile->{build},
    reference_build_source => $reference_build_profile->{source},
    remote_plot_data_policy => $keep_remote_plot_data ? 'keep_and_reuse_when_present' : 'delete_after_run_when_cleanup_enabled',
    generated_configs => [
        grep { defined $_ && length $_ } (
            $generated->{merge_config},
            $generated->{diff_config},
            $generated->{preset_config},
            $generated->{runner_config},
        )
    ],
);

if ($mode eq 'configs') {
    print_summary(\%summary);
    exit 0;
}

my %wanted = map { $_ => 1 } grep { length } split /\s*,\s*/, $plots;
my @plot_sequence = grep { $wanted{$_} } qw(manhattan local_manhattan local_gtf forest);
my $share_remote_data = (!$skip_plots && @plot_sequence > 1) ? 1 : 0;
my $shared_remote_basename = basename($runner_cfg->{DATA_GZ} // $generated->{wide_output});
my $plot_clean_oda_input = $share_remote_data ? 0 : $clean_oda_input;
my $plot_skip_upload = 0;
if (defined $env_clean_oda_input) {
    $plot_clean_oda_input = $env_clean_oda_input;
}
if (defined $env_skip_data_upload) {
    $plot_skip_upload = $env_skip_data_upload;
}
my @step_defs;

if ($source_mode eq 'raw_pgc_vcf_sumstats') {
    push @step_defs,
      {
        name        => 'merge_raw',
        description => 'Merge raw GWAS files into one long normalized table',
        command     => qq{perl "$deps_dir/merge_pgc_vcf_sumstats_long.pl" --config "$generated->{merge_config}"},
        outputs     => [ $generated->{merge_output}, $generated->{merge_manifest} ],
      },
      {
        name        => 'sort_long',
        description => 'Sort merged long GWAS by coordinate and create bgzip/tabix outputs',
        command     => qq{"$bash_path" -lc 'cd "$workdir" && INPUT_GZ="$generated->{merge_output}" OUTPUT_GZ="$generated->{sorted_output}" EXCLUDED_GZ="$generated->{excluded_output}" TMPDIR_SORT="$generated->{sort_tmpdir}" "$deps_dir/sort_long_gwas_by_coord.sh"'},
        outputs     => [ $generated->{sorted_output}, "$generated->{sorted_output}.tbi", $generated->{excluded_output} ],
      },
      {
        name        => 'diff_pairs',
        description => 'Compute pairwise differential GWAS effects',
        command     => qq{perl "$deps_dir/diff_pairwise_gwas.pl" --config "$generated->{diff_config}"},
        outputs     => [ $generated->{diff_output}, $generated->{diff_manifest} ],
      },
      {
        name        => 'standardize_diff',
        description => 'Standardize differential Z-scores and P-values',
        command     => qq{perl "$deps_dir/standardize_diff_gwas_zscore.pl" --input "$generated->{diff_output}" --output "$generated->{stdized_output}" --manifest "$generated->{stdized_manifest}" --z-col DIFF_Z},
        outputs     => [ $generated->{stdized_output}, $generated->{stdized_manifest} ],
      };
}
elsif ($source_mode eq 'precomputed_diff') {
    die "input_diff is required for source_mode=precomputed_diff\n" unless length cfg_or($spec, 'input_diff', '');
    $generated->{stdized_output} = normalize_unix_path(cfg_or($spec, 'stdized_output', $generated->{stdized_output}));
    $generated->{stdized_manifest} = normalize_unix_path(cfg_or($spec, 'stdized_manifest', $generated->{stdized_manifest}));
    push @step_defs,
      {
        name        => 'standardize_diff',
        description => 'Standardize precomputed differential GWAS results',
        command     => qq{perl "$deps_dir/standardize_diff_gwas_zscore.pl" --input "} . normalize_unix_path($spec->{input_diff}) . qq{" --output "$generated->{stdized_output}" --manifest "$generated->{stdized_manifest}" --z-col DIFF_Z},
        outputs     => [ $generated->{stdized_output}, $generated->{stdized_manifest} ],
      };
}
elsif ($source_mode eq 'precomputed_diff_stdized') {
    $generated->{stdized_output} = normalize_unix_path($spec->{input_stdized});
    die "input_stdized is required for source_mode=precomputed_diff_stdized\n"
      unless length $generated->{stdized_output};
}
elsif ($source_mode eq 'merged_gwas_table') {
    die "input_merged is required for source_mode=merged_gwas_table\n"
      unless length cfg_or($spec, 'input_merged', '');
}

if ($source_mode eq 'merged_gwas_table') {
    push @step_defs,
      {
        name        => 'extract_wide_subset',
        description => 'Normalize merged-wide GWAS results into the plotting-wide schema',
        command     => qq{perl "$deps_dir/convert_merged_gwas_to_plotwide.pl" --config "$generated->{preset_config}"},
        outputs     => [ $generated->{wide_output}, $generated->{wide_manifest} ],
      };
}
else {
    push @step_defs,
      {
        name        => 'extract_wide_subset',
        description => 'Extract the plotting-ready wide differential GWAS subset',
        command     => qq{perl "$deps_dir/extract_significant_diff_gwas.pl" --config "$generated->{preset_config}"},
        outputs     => [ $generated->{wide_output}, $generated->{wide_manifest} ],
      };
}

if (!$skip_plots) {
    if ($share_remote_data) {
        $summary{shared_remote_plot_data} = 'yes';
        $summary{shared_remote_plot_basename} = $shared_remote_basename;
        print "[info] upload_shared_plot_data deferred to plot runners for reliability.\n";
        $summary{shared_remote_plot_data_reused} = 'runner_managed_upload_or_reuse';
    } else {
        $summary{shared_remote_plot_data} = 'no';
        $summary{shared_remote_plot_data_reused} = 'not_applicable_single_plot_stage';
    }

    push @step_defs,
      {
        name        => 'plot_manhattan',
        description => 'Run the genome-wide Manhattan SAS ODA plot',
        command     => qq{"$bash_path" -lc 'cd "$workdir" && RUNNER_CONFIG_JSON="$generated->{runner_config}" SESSION_ID="$runner_session" OPEN_RESULT="$open_result" CLEAN_ODA_INPUT="$plot_clean_oda_input" SKIP_DATA_UPLOAD="$plot_skip_upload" KEEP_REMOTE_PLOT_DATA="$keep_remote_plot_data" EMIT_LOCAL_SAS_DEBUG="$emit_local_sas_scripts" LOCAL_SAS_DEBUG_ONLY="$local_sas_only" "$deps_dir/run_sas_oda_manhattan4diffgwas_download_png.sh"'},
        outputs     => $local_sas_only ? [] : [
            "$Bin/" . ($runner_cfg->{OUTPUT_PREFIX} || "${project_tag}_SAS_manhattan") . ".png",
            "$Bin/" . ($runner_cfg->{OUTPUT_PREFIX} || "${project_tag}_SAS_manhattan") . "_png.html",
        ],
        enabled     => $wanted{manhattan} ? 1 : 0,
      },
      {
        name        => 'plot_local_manhattan',
        description => 'Run the local top-hit Manhattan SAS ODA plot',
        command     => qq{"$bash_path" -lc 'cd "$workdir" && RUNNER_CONFIG_JSON="$generated->{runner_config}" SESSION_ID="$runner_session" OPEN_RESULT="$open_result" CLEAN_ODA_INPUT="$plot_clean_oda_input" SKIP_DATA_UPLOAD="$plot_skip_upload" KEEP_REMOTE_PLOT_DATA="$keep_remote_plot_data" EMIT_LOCAL_SAS_DEBUG="$emit_local_sas_scripts" LOCAL_SAS_DEBUG_ONLY="$local_sas_only" "$deps_dir/run_sas_oda_local_top_hits_manhattan_download_png.sh"'},
        outputs     => $local_sas_only ? [] : [
            "$Bin/" . ($runner_cfg->{LOCAL_OUTPUT_PREFIX} || "${project_tag}_SAS_local_top_hits_manhattan") . ".png",
            "$Bin/" . ($runner_cfg->{LOCAL_OUTPUT_PREFIX} || "${project_tag}_SAS_local_top_hits_manhattan") . ".html",
            "$Bin/" . ($runner_cfg->{LOCAL_TOP_HITS_CSV_BASENAME} || (($runner_cfg->{LOCAL_OUTPUT_PREFIX} || "${project_tag}_SAS_local_top_hits_manhattan") . "_top_hits.csv")),
        ],
        enabled     => $wanted{local_manhattan} ? 1 : 0,
      },
      {
        name        => 'plot_local_gtf',
        description => 'Run the local GTF-backed SAS ODA plot',
        command     => qq{"$bash_path" -lc 'cd "$workdir" && RUNNER_CONFIG_JSON="$generated->{runner_config}" SESSION_ID="$runner_session" OPEN_RESULT="$open_result" CLEAN_ODA_INPUT="$plot_clean_oda_input" SKIP_DATA_UPLOAD="$plot_skip_upload" KEEP_REMOTE_PLOT_DATA="$keep_remote_plot_data" EMIT_LOCAL_SAS_DEBUG="$emit_local_sas_scripts" LOCAL_SAS_DEBUG_ONLY="$local_sas_only" "$deps_dir/run_sas_oda_local_top_hits_with_gtf_download_html.sh"'},
        outputs     => $local_sas_only ? [] : [
            "$Bin/" . ($runner_cfg->{OUTPUT_HTML_BASENAME} || "${project_tag}_SAS_local_top_hits_with_gtf.html"),
            "$Bin/" . ($runner_cfg->{LOCAL_TOP_HITS_CSV_BASENAME} || (($runner_cfg->{LOCAL_OUTPUT_PREFIX} || "${project_tag}_SAS_local_top_hits_manhattan") . "_top_hits.csv")),
        ],
        enabled     => $wanted{local_gtf} ? 1 : 0,
      },
      {
        name        => 'plot_forest',
        description => 'Run the top-hit forest SAS ODA plot',
        command     => qq{"$bash_path" -lc 'cd "$workdir" && RUNNER_CONFIG_JSON="$generated->{runner_config}" SESSION_ID="$runner_session" OPEN_RESULT="$open_result" CLEAN_ODA_INPUT="$plot_clean_oda_input" SKIP_DATA_UPLOAD="$plot_skip_upload" KEEP_REMOTE_PLOT_DATA="$keep_remote_plot_data" EMIT_LOCAL_SAS_DEBUG="$emit_local_sas_scripts" LOCAL_SAS_DEBUG_ONLY="$local_sas_only" "$deps_dir/run_sas_oda_top_hits_forest_plot_download_html.sh"'},
        outputs     => $local_sas_only ? [] : [
            "$Bin/" . ($runner_cfg->{FOREST_OUTPUT_HTML_BASENAME} || (($runner_cfg->{FOREST_OUTPUT_PREFIX} || "${project_tag}_SAS_top_hits_forest") . ".html")),
            "$Bin/" . ($runner_cfg->{FOREST_TOP_HITS_CSV_BASENAME} || (($runner_cfg->{FOREST_OUTPUT_PREFIX} || "${project_tag}_SAS_top_hits_forest") . "_top_hits.csv")),
            "$Bin/" . ($runner_cfg->{FOREST_OUTPUT_MANIFEST_BASENAME} || (($runner_cfg->{FOREST_OUTPUT_PREFIX} || "${project_tag}_SAS_top_hits_forest") . ".manifest.tsv")),
        ],
        enabled     => $wanted{forest} ? 1 : 0,
      },
      {
        name        => 'cleanup_shared_plot_data',
        description => 'Delete the shared wide plot subset from SAS ODA after multi-plot runs',
        command     => qq{"$bash_path" -lc 'cd "$workdir" && perl "$oda_helper_unix" --delete-file "$shared_remote_basename" --persistent --session-id "$runner_session" --output-prefix "cleanup_shared_plot_data_} . strftime('%Y%m%d_%H%M%S', localtime) . qq{"'},
        force       => 1,
        enabled     => ($share_remote_data && $clean_oda_input && !$keep_remote_plot_data) ? 1 : 0,
      };
}

my $step_selection = resolve_step_selection(
    step_args     => \@step_args,
    step_flag     => \%step_flag,
    from_step     => $from_step,
    to_step       => $to_step,
    step_defs     => \@step_defs,
);

if ($list_steps) {
    print_step_catalog(\@step_defs, $step_selection);
    exit 0;
}

$summary{steps_selected} = [ @{ $step_selection->{selected_order} } ];
$summary{step_selection_mode} = $step_selection->{selection_mode};

for my $step (@step_defs) {
    next unless $step_selection->{selected}{$step->{name}};
    if (exists $step->{enabled} && !$step->{enabled}
        && !($step_selection->{explicit_requested}{$step->{name}})) {
        print "[skip] $step->{name} is not active for the current plots/source_mode selection\n";
        next;
    }
    if ($step->{name} =~ /^plot_/) {
        print STDERR "[info] Validating generated files for plotting step '$step->{name}'...\n";
        validate_generated_files($generated, $pair_info);
    }
    #Only skip the step when all output were found;
    foreach my $of (@{$step->{outputs}}) {
        $force = 1 unless -f $of;
    }
    run_step(
        name    => $step->{name},
        command => $step->{command},
        outputs => $step->{outputs} || [],
        force   => exists $step->{force} ? $step->{force} : $force,
    );
}

$summary{wide_output} = $generated->{wide_output};
$summary{runner_config} = $generated->{runner_config};
$summary{plots_requested} = $skip_plots ? 'none' : $plots;

# If common-association mode requested, run an independent verifier to
# enumerate all candidate common-association loci and write an easy-to-read
# TSV for inspection and downstream use. This helps catch cases where the
# SAS top-hit selection may have been overly restrictive.
if ($runner_cfg->{TOP_HIT_MODE} && lc($runner_cfg->{TOP_HIT_MODE}) eq 'common_association') {
    my $verify = File::Spec->catfile($deps_dir, 'verify_common_association_loci.pl');
    if (-e cygpath_to_win($verify)) {
        my $out_base = safe_name($artifact_stem) . '.common_assoc_verify';
        my $out_tsv = File::Spec->catfile($output_dir, $out_base . '.tsv');
        my $cand_tsv = File::Spec->catfile($output_dir, $out_base . '.candidates.tsv');
        my $thr_arg = defined $runner_cfg->{TOP_HIT_SIGNAL_THRSHDS} ? "--top-p-thresholds '$runner_cfg->{TOP_HIT_SIGNAL_THRSHDS}'" : '';
        my $cmd = qq{perl '$verify' --spec '$spec_file' --output '$out_tsv' --candidates-out '$cand_tsv' $thr_arg};
        print "[info] Running common-association verifier: $cmd\n";
        system($cmd) == 0 or warn "verify_common_association_loci.pl failed: $?\n";
        $summary{common_association_verify} = $out_tsv if -s cygpath_to_win($out_tsv);
    }
    else {
        warn "Verifier not found: $verify\n";
    }
}

print_summary(\%summary);

sub build_generated_paths {
    my (%args) = @_;
    my $configs_dir = $args{configs_dir};
    my $output_dir = $args{output_dir};
    my $stem = safe_name($args{artifact_stem});
    return {
        merge_config     => "$configs_dir/auto_${stem}_merge.json",
        diff_config      => "$configs_dir/auto_${stem}_diff.json",
        preset_config    => "$configs_dir/auto_${stem}_preset.json",
        runner_config    => "$configs_dir/auto_${stem}_runner.json",
        merge_output     => "$output_dir/${stem}_merged_long.tsv.gz",
        merge_manifest   => "$output_dir/${stem}_merged_long.manifest.tsv",
        sorted_output    => "$output_dir/${stem}_merged_long.sorted.coord.tsv.gz",
        excluded_output  => "$output_dir/${stem}_merged_long.sorted.excluded_noncoord.tsv.gz",
        sort_tmpdir      => "$output_dir/sort_tmp_${stem}",
        diff_output      => "$output_dir/${stem}.tsv.gz",
        diff_manifest    => "$output_dir/${stem}.manifest.tsv",
        stdized_output   => "$output_dir/${stem}.stdized.tsv.gz",
        stdized_manifest => "$output_dir/${stem}.stdized.manifest.tsv",
        wide_output      => "$output_dir/${stem}.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz",
        wide_manifest    => "$output_dir/${stem}.stdized.wide_beta_se_p_p_lt_0p05.final.manifest.tsv",
    };
}

sub resolve_reference_build_profile_for_spec {
    my (%args) = @_;
    my $spec = $args{spec} || {};
    my $explicit = $args{reference_override};
    $explicit = cfg_or($spec, 'reference_build',
        cfg_or($spec, 'genome_build',
            cfg_or($spec, 'gtf_reference_build', '')))
      unless defined $explicit && length $explicit;

    my @paths = collect_reference_build_candidate_files(
        spec        => $spec,
        source_mode => $args{source_mode},
    );
    my @headers;
    for my $path (@paths) {
        my $header = read_first_line_maybe_gzip($path);
        push @headers, $header if defined $header && length $header;
        last if @headers >= 3;
    }

    my $profile = detect_reference_build_profile(
        explicit_build => $explicit,
        file_paths     => \@paths,
        header_lines   => \@headers,
        default_build  => 'hg38',
    );

    if (($profile->{source} || '') eq 'fallback_default') {
        print STDERR "WARNING: Reference build could not be confidently inferred from the spec inputs; defaulting to hg38/GRCh38. "
          . "Set reference_build in the spec or pass --reference-build hg19|hg38|t2t to override.\n";
    }
    else {
        print "[info] Resolved reference build: $profile->{build} ($profile->{source}";
        print ", $profile->{evidence}" if defined $profile->{evidence} && length $profile->{evidence};
        print ")\n";
    }

    return $profile;
}

sub collect_reference_build_candidate_files {
    my (%args) = @_;
    my $spec = $args{spec} || {};
    my $input_dir = normalize_unix_path(cfg_or($spec, 'input_dir', ''));
    my @paths;
    my %seen;

    my $push_path = sub {
        my ($path) = @_;
        return unless defined $path && length $path;
        $path = normalize_unix_path($path);
        if ($path !~ m{^(?:[A-Za-z]:)?[\\/]} && $input_dir) {
            $path = normalize_unix_path(File::Spec->catfile($input_dir, $path));
        }
        return if $seen{$path}++;
        push @paths, $path;
    };

    for my $key (qw(
        source_long_gz source_long long_gz long_tsv input_gz input_tsv
        precomputed_long_gz precomputed_long_tsv precomputed_wide_gz precomputed_wide_tsv
        wide_gz wide_tsv data_gz input_merged
    )) {
        next unless exists $spec->{$key};
        if (ref($spec->{$key}) eq 'ARRAY') {
            $push_path->($_) for @{ $spec->{$key} || [] };
        }
        else {
            $push_path->($spec->{$key});
        }
    }

    for my $group (@{ $spec->{groups} || [] }) {
        for my $file (@{ $group->{files} || [] }) {
            $push_path->($file);
        }
    }

    return grep { defined $_ && length $_ } @paths;
}

sub read_first_line_maybe_gzip {
    my ($path) = @_;
    return '' unless defined $path && length $path && -f $path;

    if ($path =~ /\.gz$/i) {
        my $fh = IO::Uncompress::Gunzip->new($path);
        return '' unless $fh;
        my $line = <$fh>;
        close $fh;
        return sanitize_header_line($line);
    }

    open my $fh, '<', $path or return '';
    my $line = <$fh>;
    close $fh;
    return sanitize_header_line($line);
}

sub sanitize_header_line {
    my ($line) = @_;
    return '' unless defined $line;
    chomp $line;
    $line =~ s/\r$//;
    return $line;
}

sub infer_spec_from_gwas_dir {
    my (%args) = @_;
    my $gwas_dir = normalize_unix_path($args{gwas_dir} // '');
    die "--gwas-dir is required for auto spec generation\n" unless length $gwas_dir;
    my $workdir = normalize_unix_path($args{workdir} // script_root_dir());
    my $configs_dir = "$workdir/configs";

    progress_note("Scanning GWAS directory: $gwas_dir");
    my $study_prefix = infer_study_prefix_from_dir($gwas_dir);
    my @raw_files = discover_raw_gwas_candidate_files($gwas_dir);
    my @precomputed_files = discover_precomputed_candidate_files($gwas_dir);
    my @merged_files = discover_merged_candidate_files($gwas_dir);
    progress_note(
        "Found candidate files: raw=" . scalar(@raw_files)
        . ", precomputed=" . scalar(@precomputed_files)
        . ", merged=" . scalar(@merged_files)
    );
    die "No candidate GWAS summary-statistics files found in $gwas_dir\n"
      unless @raw_files || @precomputed_files || @merged_files;

    if (@precomputed_files) {
        progress_note("Checking precomputed differential GWAS candidates...");
        my ($pre_spec, $pre_path) = infer_precomputed_spec_from_dir(
            gwas_dir                => $gwas_dir,
            workdir                 => $workdir,
            configs_dir             => $configs_dir,
            study_prefix            => $study_prefix,
            project_tag_override    => $args{project_tag_override},
            artifact_stem_override  => $args{artifact_stem_override},
            spec_out                => $args{spec_out},
            files                   => \@precomputed_files,
        );
        return ($pre_spec, $pre_path) if $pre_spec;
    }

    if (@merged_files) {
        progress_note("Checking merged wide GWAS candidates...");
        my ($merged_spec, $merged_path) = infer_merged_spec_from_dir(
            gwas_dir                => $gwas_dir,
            workdir                 => $workdir,
            configs_dir             => $configs_dir,
            study_prefix            => $study_prefix,
            project_tag_override    => $args{project_tag_override},
            artifact_stem_override  => $args{artifact_stem_override},
            spec_out                => $args{spec_out},
            files                   => \@merged_files,
        );
        return ($merged_spec, $merged_path) if $merged_spec;
    }

    my (%groups_by_tag, @detected_files, @skipped_files);
    my (%effect_metric_seen, %format_class_seen);
    my $raw_total = scalar @raw_files;
    my $raw_idx = 0;
    for my $path (@raw_files) {
        $raw_idx++;
        progress_note("Inspecting raw GWAS [$raw_idx/$raw_total]: " . basename($path));
        my $info = eval { inspect_raw_gwas_file($path, $args{raw_column_aliases}) };
        if ($@) {
            chomp(my $err = $@);
            push @skipped_files, { file => basename($path), reason => $err };
            progress_note("Skipped " . basename($path) . ": $err");
            next;
        }

        my $tag = infer_group_tag_from_file($path, $study_prefix);
        my $label = infer_group_label_from_file($path, $tag);
        $groups_by_tag{$tag} ||= {
            tag   => $tag,
            label => $label,
            files => [],
            format_class  => $info->{format_class},
            effect_metric => $info->{effect_metric},
        };
        push @{ $groups_by_tag{$tag}{files} }, basename($path);
        push @detected_files, {
            file             => basename($path),
            inferred_group   => $tag,
            inferred_label   => $label,
            effect_metric    => $info->{effect_metric},
            format_class     => $info->{format_class},
            resolved_columns => $info->{resolved_columns},
        };
        $effect_metric_seen{ $info->{effect_metric} || 'UNKNOWN' }++;
        $format_class_seen{ $info->{format_class} || 'UNKNOWN' }++;
        progress_note(
            "Accepted " . basename($path)
            . " as $tag"
            . " [" . ($info->{format_class} || 'UNKNOWN')
            . "; " . ($info->{effect_metric} || 'UNKNOWN') . "]"
        );
    }

    my @all_groups = sort { $a->{tag} cmp $b->{tag} } values %groups_by_tag;
    my @selection_warnings;
    my @groups = select_compatible_groups(\@all_groups, \@selection_warnings);
    progress_note(
        "Raw GWAS grouping complete: accepted_groups=" . scalar(@all_groups)
        . ", selected_groups=" . scalar(@groups)
        . ", skipped_files=" . scalar(@skipped_files)
    );
    if (@groups < 2) {
        progress_note("Not enough compatible raw groups; retrying precomputed diff detection...");
        my ($pre_spec, $pre_path) = infer_precomputed_spec_from_dir(
            gwas_dir                => $gwas_dir,
            workdir                 => $workdir,
            configs_dir             => $configs_dir,
            study_prefix            => $study_prefix,
            project_tag_override    => $args{project_tag_override},
            artifact_stem_override  => $args{artifact_stem_override},
            spec_out                => $args{spec_out},
            files                   => \@precomputed_files,
        );
        return ($pre_spec, $pre_path) if $pre_spec;
        die "Need at least two valid raw GWAS groups or one valid precomputed diff/stdized table in $gwas_dir\n";
    }

    my $project_tag = safe_name($args{project_tag_override} || infer_project_tag_from_dir($gwas_dir, \@groups));
    my $artifact_stem = safe_name($args{artifact_stem_override} || ($project_tag . '_effects'));

    my @pairs;
    my %prefix_seen;
    for my $i (0 .. $#groups - 1) {
        for my $j ($i + 1 .. $#groups) {
            my $g1 = $groups[$i];
            my $g2 = $groups[$j];
            my $prefix = safe_name($g1->{label}) . '_' . safe_name($g2->{label});
            $prefix = uc $prefix;
            my $base_prefix = $prefix;
            my $suffix = 2;
            while ($prefix_seen{$prefix}) {
                $prefix = $base_prefix . '_' . $suffix++;
            }
            $prefix_seen{$prefix} = 1;
            push @pairs, {
                pair_tag      => $g1->{tag} . '_vs_' . $g2->{tag},
                group1        => $g1->{tag},
                group2        => $g2->{tag},
                group1_label  => $g1->{label},
                group2_label  => $g2->{label},
                prefix        => $prefix,
                label         => $g1->{label} . ' vs ' . $g2->{label},
            };
        }
    }

    my $focus_prefix = $pairs[0]{prefix};
    my @warnings;
    if (keys(%effect_metric_seen) > 1) {
        push @warnings, 'Mixed effect metrics detected across raw GWAS files: '
          . join(', ', sort keys %effect_metric_seen)
          . '. The merge step currently converts OR/odds-ratio inputs to log(OR), but comparisons across mixed source formats should be reviewed carefully.';
    }
    if (keys(%format_class_seen) > 1) {
        push @warnings, 'Mixed raw GWAS schema families detected: '
          . join(', ', sort keys %format_class_seen)
          . '. Auto-generated comparisons may be technically runnable but should be checked for scientific comparability.';
    }
    push @warnings, @selection_warnings if @selection_warnings;
    my $detected_build = detect_reference_build_profile(
        file_paths    => \@raw_files,
        default_build => 'hg38',
    );
    my $spec = {
        source_mode            => 'raw_pgc_vcf_sumstats',
        project_tag            => $project_tag,
        artifact_stem          => $artifact_stem,
        reference_build        => ($detected_build->{build} || 'hg38'),
        input_dir              => $gwas_dir,
        output_dir             => $gwas_dir,
        workdir                => $workdir,
        cygwin_bash            => '/bin/bash',
        threshold              => 0.05,
        keep_remote_plot_data  => 1,
        top_hit_focus_prefix   => $focus_prefix,
        top_hit_signal_thrshd  => '1e-6',
        top_hit_dist_bp        => '1e6',
        local_window_bp        => '1e7',
        local_gtf_window_bp    => '1e7',
        include_non_protein_coding_genes_in_local_gtf => 0,
        open_result            => 0,
        clean_oda_input        => 1,
        (ref($args{raw_column_aliases}) eq 'HASH' && keys %{ $args{raw_column_aliases} || {} }
            ? (raw_column_aliases => $args{raw_column_aliases})
            : ()),
        groups                 => [ map { { tag => $_->{tag}, files => $_->{files} } } @groups ],
        pairs                  => \@pairs,
        auto_detect_summary    => {
            generated_from_dir => $gwas_dir,
            detected_files     => \@detected_files,
            skipped_files      => \@skipped_files,
            warnings           => \@warnings,
        },
    };

    my $spec_path = normalize_unix_path($args{spec_out} || "$configs_dir/auto_${artifact_stem}_from_dir.spec.json");
    progress_note(
        "Draft spec ready: groups=" . scalar(@groups)
        . ", pairs=" . scalar(@pairs)
        . ", output=$spec_path"
    );
    return ($spec, $spec_path);
}

sub discover_raw_gwas_candidate_files {
    my ($dir) = @_;
    my $win_dir = cygpath_to_win($dir);
    die "GWAS directory does not exist: $dir\n" unless -d $win_dir;

    opendir(my $dh, $win_dir) or die "Cannot open GWAS directory $win_dir: $!\n";
    my @files;
    while (my $entry = readdir($dh)) {
        next if $entry eq '.' || $entry eq '..';
        my $full = "$win_dir\\$entry";
        next unless -f $full;
        next if $entry =~ /(?:merged_long|diff_effects|\.stdized\b|wide_beta_se_p|excluded_noncoord|sorted\.coord|manifest|merged_plotwide|gunplot|SAS_top_hits_forest|single_snp\.data|\.gp$|\.html$|\.png$)/i;
        next unless $entry =~ /\.(?:tsv|txt|sumstats|assoc|vcf)(?:\.(?:gz|bgz|bgzip))?$/i
                 || $entry =~ /\.gz$/i
                 || $entry =~ /\.vcf\.tsv\.gz$/i;
        push @files, win_to_cygpath($full);
    }
    closedir $dh;
    return sort @files;
}

sub discover_precomputed_candidate_files {
    my ($dir) = @_;
    my $win_dir = cygpath_to_win($dir);
    die "GWAS directory does not exist: $dir\n" unless -d $win_dir;

    opendir(my $dh, $win_dir) or die "Cannot open GWAS directory $win_dir: $!\n";
    my @files;
    while (my $entry = readdir($dh)) {
        next if $entry eq '.' || $entry eq '..';
        my $full = "$win_dir\\$entry";
        next unless -f $full;
        next if $entry =~ /\.vcf\.tsv\.gz$/i;
        next unless $entry =~ /\.(?:tsv|txt)(?:\.(?:gz|bgz|bgzip))?$/i;
        next if $entry =~ /(?:manifest|excluded_noncoord|sorted\.coord\.tsv\.gz\.tbi|merged_plotwide|gunplot|SAS_top_hits_forest|single_snp\.data|\.gp$|\.html$|\.png$)/i;
        push @files, win_to_cygpath($full);
    }
    closedir $dh;
    return sort {
        precomputed_priority($a) <=> precomputed_priority($b)
          || $a cmp $b
    } @files;
}

sub discover_merged_candidate_files {
    my ($dir) = @_;
    my $win_dir = cygpath_to_win($dir);
    die "GWAS directory does not exist: $dir\n" unless -d $win_dir;

    opendir(my $dh, $win_dir) or die "Cannot open GWAS directory $win_dir: $!\n";
    my @files;
    while (my $entry = readdir($dh)) {
        next if $entry eq '.' || $entry eq '..';
        my $full = "$win_dir\\$entry";
        next unless -f $full;
        next unless $entry =~ /(?:merged|plus_meta|meta|combined|wide)/i;
        next unless $entry =~ /\.(?:tsv|txt)(?:\.(?:gz|bgz|bgzip))?$/i
                 || $entry =~ /\.gz$/i;
        next if $entry =~ /(?:manifest|merged_plotwide|gunplot|SAS_top_hits_forest|single_snp\.data|\.gp$|\.html$|\.png$)/i;
        push @files, win_to_cygpath($full);
    }
    closedir $dh;
    return sort @files;
}

sub precomputed_priority {
    my ($path) = @_;
    my $base = lc basename($path);
    return 0 if $base =~ /\.stdized\.tsv\.gz$/;
    return 1 if $base =~ /diff_effects.*\.tsv\.gz$/;
    return 2 if $base =~ /\.tsv\.gz$/;
    return 3;
}

sub inspect_merged_wide_gwas_file {
    my ($path) = @_;
    my $fh = open_text_reader($path);

    my ($header, @cols);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^##/;
        next if $line =~ /^\s*$/;
        $header = $line;
        @cols = split /\t/, $header, -1;
        last;
    }
    close $fh;
    die "No readable header found in $path\n" unless defined $header;

    my %colset = map { $_ => 1 } @cols;
    for my $required (qw(CHR BP SNP)) {
        my $actual = find_header_actual(\@cols, [$required]);
        die "Required merged-wide column $required not found in $path\n" unless defined $actual;
    }
    my $a1_col = find_header_actual(\@cols, [qw(A1 EA EFFECT_ALLELE ALT)]);
    my $a2_col = find_header_actual(\@cols, [qw(A2 NEA OTHER_ALLELE REF)]);

    my @ordered_group_stems;
    my %group_track;
    for my $col (@cols) {
        next unless $col =~ /^BETA_(.+)$/i;
        my $stem = $1;
        my $se_col = find_header_actual(\@cols, ["SE_$stem"]);
        my $p_col = find_header_actual(\@cols, ["P_$stem"]);
        next unless defined $se_col && defined $p_col;
        my $norm_stem = normalize_merged_track_id($stem);
        next unless length $norm_stem;
        next if exists $group_track{$norm_stem};
        $group_track{$norm_stem} = {
            id       => $norm_stem,
            label    => $stem,
            beta_col => find_header_actual(\@cols, ["BETA_$stem"]),
            se_col   => $se_col,
            p_col    => $p_col,
        };
        push @ordered_group_stems, $norm_stem;
    }
    for my $col (@cols) {
        next unless $col =~ /^(.+)_BETA$/i;
        my $stem = $1;
        my $se_col = find_header_actual(\@cols, ["${stem}_SE"]);
        my $p_col = find_header_actual(\@cols, ["${stem}_P"]);
        next unless defined $se_col && defined $p_col;
        my $norm_stem = normalize_merged_track_id($stem);
        next unless length $norm_stem;
        next if exists $group_track{$norm_stem};
        $group_track{$norm_stem} = {
            id       => $norm_stem,
            label    => $stem,
            beta_col => find_header_actual(\@cols, ["${stem}_BETA"]),
            se_col   => $se_col,
            p_col    => $p_col,
        };
        push @ordered_group_stems, $norm_stem;
    }

    my @group_tracks = map { $group_track{$_} } @ordered_group_stems;
    die "Could not identify at least two cohort-level BETA/SE/P track blocks in $path\n"
      unless @group_tracks >= 2;

    my %group_source_cols = map {
        $_->{beta_col} => 1,
        $_->{se_col}   => 1,
        $_->{p_col}    => 1,
    } @group_tracks;
    my @extra_tracks;
    my %seen_extra;
    for my $col (@cols) {
        next if $group_source_cols{$col};
        my ($stem, $z_col, $p_col);
        if ($col =~ /^WEIGHTED_Z_(.+)$/i) {
            $stem = $1;
            $z_col = find_header_actual(\@cols, ["WEIGHTED_Z_$stem"]);
            $p_col = find_header_actual(\@cols, ["PR_$stem", "PWZ_$stem", "P_$stem", "${stem}_P"]);
        }
        elsif ($col =~ /^(.+)_Z$/i) {
            $stem = $1;
            $z_col = find_header_actual(\@cols, ["${stem}_Z"]);
            $p_col = find_header_actual(\@cols, ["${stem}_P", "P_$stem"]);
        }
        next unless defined $stem && defined $z_col && defined $p_col;
        my $norm_id = normalize_merged_track_id($stem);
        next unless length $norm_id;
        next if $group_track{$norm_id};
        next if $seen_extra{$norm_id}++;
        push @extra_tracks, {
            id       => $norm_id,
            label    => $stem,
            p_col    => $p_col,
            z_col    => $z_col,
            beta_col => find_header_actual(\@cols, ["BETA_$stem", "${stem}_BETA", "BETAR_$stem", "${stem}_BETA"]),
            se_col   => find_header_actual(\@cols, ["SE_$stem", "${stem}_SE"]),
        };
    }

    return {
        header       => $header,
        cols         => \@cols,
        chr_col      => find_header_actual(\@cols, [qw(CHR chromosome chrom)]),
        bp_col       => find_header_actual(\@cols, [qw(BP POS position)]),
        snp_col      => find_header_actual(\@cols, [qw(SNP rsid marker id)]),
        a1_col       => $a1_col,
        a2_col       => $a2_col,
        group_tracks => \@group_tracks,
        extra_tracks => \@extra_tracks,
    };
}

sub infer_merged_spec_from_dir {
    my (%args) = @_;
    my @files = @{ $args{files} || [] };
    my @candidates;
    my $total = scalar @files;
    my $idx = 0;
    for my $path (@files) {
        $idx++;
        progress_note("Inspecting merged-wide GWAS [$idx/$total]: " . basename($path));
        my $info = eval { inspect_merged_wide_gwas_file($path) };
        if ($@) {
            chomp(my $err = $@);
            progress_note("Skipped merged-wide candidate " . basename($path) . ": $err");
            next;
        }
        push @candidates, { path => $path, %{$info} };
        progress_note(
            "Accepted merged-wide candidate " . basename($path)
            . " with cohort_tracks=" . scalar(@{ $info->{group_tracks} || [] })
            . " and extra_tracks=" . scalar(@{ $info->{extra_tracks} || [] })
        );
    }
    return undef unless @candidates;

    @candidates = sort {
        scalar(@{ $b->{group_tracks} || [] }) <=> scalar(@{ $a->{group_tracks} || [] })
        || scalar(@{ $b->{extra_tracks} || [] }) <=> scalar(@{ $a->{extra_tracks} || [] })
        || length($b->{path}) <=> length($a->{path})
    } @candidates;
    my $best = $candidates[0];

    my @pairs = build_pairs_from_group_track_labels(
        $best->{group_tracks},
        $args{study_prefix},
    );
    die "Could not infer at least one cohort comparison from merged wide file $best->{path}\n"
      unless @pairs;

    my $project_tag = safe_name($args{project_tag_override} || infer_project_tag_from_dir($args{gwas_dir}, []));
    my $artifact_stem = safe_name($args{artifact_stem_override} || ($project_tag . '_merged'));
    my $detected_build = detect_reference_build_profile(
        file_paths    => [ $best->{path} ],
        header_lines  => [ $best->{header} ],
        default_build => 'hg38',
    );
    my $spec = {
        source_mode           => 'merged_gwas_table',
        project_tag           => $project_tag,
        artifact_stem         => $artifact_stem,
        reference_build       => ($detected_build->{build} || 'hg38'),
        input_dir             => $args{gwas_dir},
        output_dir            => $args{gwas_dir},
        workdir               => $args{workdir},
        cygwin_bash           => '/bin/bash',
        threshold             => 0.05,
        keep_remote_plot_data => 1,
        top_hit_focus_prefix  => $pairs[0]{prefix},
        top_hit_signal_thrshd => '1e-6',
        top_hit_dist_bp       => '1e6',
        local_window_bp       => '1e7',
        local_gtf_window_bp   => '1e7',
        include_non_protein_coding_genes_in_local_gtf => 0,
        open_result           => 0,
        clean_oda_input       => 1,
        input_merged          => $best->{path},
        merged_base_cols      => {
            chr => $best->{chr_col},
            bp  => $best->{bp_col},
            snp => $best->{snp_col},
            a1  => $best->{a1_col},
            a2  => $best->{a2_col},
        },
        merged_group_tracks   => $best->{group_tracks},
        merged_extra_tracks   => $best->{extra_tracks},
        extra_tracks          => [
            map {
                +{
                    id      => $_->{id},
                    label   => normalize_merged_track_label($_->{label}),
                    kind    => 'extra',
                    pvar    => $_->{id} . '_P',
                    zvar    => $_->{id} . '_Z',
                    betavar => (defined($_->{beta_col}) && length($_->{beta_col}) ? $_->{id} . '_BETA' : ''),
                    sevar   => (defined($_->{se_col}) && length($_->{se_col}) ? $_->{id} . '_SE' : ''),
                }
            } @{ $best->{extra_tracks} || [] }
        ],
        pairs                 => \@pairs,
        auto_detect_summary   => {
            generated_from_dir => $args{gwas_dir},
            selected_file      => basename($best->{path}),
            selected_mode      => 'merged_gwas_table',
            cohort_tracks      => [ map { $_->{label} } @{ $best->{group_tracks} || [] } ],
            extra_tracks       => [ map { $_->{label} } @{ $best->{extra_tracks} || [] } ],
        },
    };
    my $spec_path = normalize_unix_path($args{spec_out} || "$args{configs_dir}/auto_${artifact_stem}_from_dir.spec.json");
    progress_note(
        "Selected merged-wide GWAS source: " . basename($best->{path})
        . "; pairs=" . scalar(@pairs)
        . ", output=$spec_path"
    );
    return ($spec, $spec_path);
}

sub build_pairs_from_group_track_labels {
    my ($tracks, $study_prefix) = @_;
    my @tracks = @{ $tracks || [] };
    my @pairs;
    my %prefix_seen;
    for my $i (0 .. $#tracks - 1) {
        for my $j ($i + 1 .. $#tracks) {
            my $g1 = $tracks[$i];
            my $g2 = $tracks[$j];
            next unless ref($g1) eq 'HASH' && ref($g2) eq 'HASH';
            my $group1_code = normalize_merged_track_id($g1->{label});
            my $group2_code = normalize_merged_track_id($g2->{label});
            next unless length($group1_code) && length($group2_code);
            my $group1_tag = compose_group_tag($study_prefix, $group1_code);
            my $group2_tag = compose_group_tag($study_prefix, $group2_code);
            my $group1_label = normalize_merged_track_label($g1->{label});
            my $group2_label = normalize_merged_track_label($g2->{label});
            my $prefix = uc(safe_name($group1_label) . '_' . safe_name($group2_label));
            my $base = $prefix;
            my $suffix = 2;
            while ($prefix_seen{$prefix}) {
                $prefix = $base . '_' . $suffix++;
            }
            $prefix_seen{$prefix} = 1;
            push @pairs, {
                pair_tag      => $group1_tag . '_vs_' . $group2_tag,
                group1        => $group1_tag,
                group2        => $group2_tag,
                group1_label  => $group1_label,
                group2_label  => $group2_label,
                prefix        => $prefix,
                label         => $group1_label . ' vs ' . $group2_label,
                source_group1 => $g1->{id},
                source_group2 => $g2->{id},
            };
        }
    }
    return @pairs;
}

sub normalize_merged_track_id {
    my ($value) = @_;
    $value = '' unless defined $value;
    $value =~ s/^\s+|\s+$//g;
    $value =~ s/[^A-Za-z0-9]+/_/g;
    $value =~ s/^_+|_+$//g;
    return uc($value);
}

sub normalize_merged_track_label {
    my ($value) = @_;
    $value = '' unless defined $value;
    $value =~ s/^\s+|\s+$//g;
    $value =~ s/\s+/ /g;
    $value =~ s/^meta$/Meta/i;
    return $value;
}

sub inspect_raw_gwas_file {
    my ($path, $alias_overrides) = @_;
    my $fh = open_text_reader($path);

    my ($header, @cols);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^##/;
        next if $line =~ /^\s*$/;
        $header = $line;
        @cols = split /\t/, $header, -1;
        last;
    }
    close $fh;

    die "No readable header found in $path\n" unless defined $header;
    my $data_rows = count_data_rows($path, 1);
    die "No data rows found after header in $path\n" unless $data_rows > 0;
    my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
    my %resolved = resolve_header_aliases_local(\@cols, \%idx, $path, $header, $alias_overrides);
    my @notes = map {
        my $actual = exists $resolved{$_} ? $cols[ $resolved{$_} ] : 'NA';
        $_ . '=' . $actual;
    } qw(CHROM ID POS A1 A2 BETA SE PVAL FCAS FCON IMPINFO NCAS NCON NEFF);
    my $beta_actual = $cols[ $resolved{BETA} ] // '';
    my $beta_norm = normalize_header_local($beta_actual);
    my $effect_metric = ($beta_norm eq 'OR' || $beta_norm eq 'ODDSRATIO') ? 'LOG_OR_FROM_OR' : 'BETA_LIKE';
    my $format_class = infer_raw_format_class(\@cols, $path);

    return {
        header           => $header,
        effect_metric    => $effect_metric,
        format_class     => $format_class,
        data_rows        => $data_rows,
        resolved_columns => \@notes,
    };
}

sub count_data_rows {
    my ($path, $limit) = @_;
    $limit ||= 1;
    my $fh = open_text_reader($path);
    my $header_seen = 0;
    my $rows = 0;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^##/;
        next if $line =~ /^\s*$/;
        if (!$header_seen) {
            $header_seen = 1;
            next;
        }
        $rows++;
        last if $rows >= $limit;
    }
    close $fh;
    return $rows;
}

sub infer_raw_format_class {
    return DiffGWASRawSchema::infer_raw_format_class(@_);
}

sub inspect_precomputed_diff_file {
    my ($path) = @_;
    my $fh = open_text_reader($path);

    my ($header, @cols);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^##/;
        next if $line =~ /^\s*$/;
        if (!defined $header) {
            $header = $line;
            @cols = split /\t/, $header, -1;
            last;
        }
    }
    close $fh;
    die "No readable header found in $path\n" unless defined $header;
    my $pair_col = find_header_actual(\@cols, [qw(PAIR_TAG PAIRTAG COMPARISON COMPARE_TAG)]);
    die "PAIR_TAG column not found in $path\n" unless defined $pair_col;

    my @sample_pair_tags;
    my %seen_pair;
    $fh = open_text_reader($path);
    my $header_seen = 0;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^##/;
        next if $line =~ /^\s*$/;
        if (!$header_seen) {
            $header_seen = 1;
            next;
        }
        my @v = split /\t/, $line, -1;
        my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
        my $tag = $v[$idx{$pair_col}] // '';
        if (length $tag && !$seen_pair{$tag}++) {
            push @sample_pair_tags, $tag;
            last if @sample_pair_tags >= 50;
        }
    }
    close $fh;
    my @required = qw(CHR BP A1 A2 SNP PAIR_TAG DIFF_BETA DIFF_SE DIFF_P);
    my @std_extra = qw(STD_DIFF_Z STD_DIFF_P);

    for my $need (@required) {
        my $actual = find_header_actual(\@cols, [$need]);
        die "Required column $need not found in $path\n" unless defined $actual;
    }
    my @side_stems = detect_precomputed_side_stems(\@cols);
    die "Could not identify the two compared GWAS stems in $path\n" unless @side_stems == 2;
    for my $stem (@side_stems) {
        for my $suffix (qw(BETA SE P)) {
            my $actual = find_header_actual(\@cols, [$stem . '_' . $suffix]);
            die "Required column ${stem}_${suffix} not found in $path\n" unless defined $actual;
        }
    }
    my $has_std = 1;
    for my $need (@std_extra) {
        my $actual = find_header_actual(\@cols, [$need]);
        $has_std = 0 unless defined $actual;
    }
    my $has_diff_z = find_header_actual(\@cols, [qw(DIFF_Z DIFFZ)]) ? 1 : 0;
    my $mode = $has_std ? 'precomputed_diff_stdized' : ($has_diff_z ? 'precomputed_diff' : '');
    die "Could not classify precomputed diff file $path as diff or stdized\n" unless $mode;

    return {
        mode           => $mode,
        header         => $header,
        sample_pairs   => \@sample_pair_tags,
        side_stems     => \@side_stems,
    };
}

sub detect_precomputed_side_stems {
    my ($cols) = @_;
    my %stem_ok;
    for my $col (@{$cols}) {
        next unless $col =~ /^(.+)_BETA$/;
        my $stem = $1;
        next if $stem =~ /^(?:DIFF|STD_DIFF|ORIG_DIFF)$/;
        $stem_ok{$stem}{BETA} = 1;
    }
    for my $col (@{$cols}) {
        for my $suffix (qw(SE P)) {
            if ($col =~ /^(.+)_\Q$suffix\E$/) {
                my $stem = $1;
                $stem_ok{$stem}{$suffix} = 1 if exists $stem_ok{$stem};
            }
        }
    }
    my @stems = grep { $stem_ok{$_}{BETA} && $stem_ok{$_}{SE} && $stem_ok{$_}{P} } sort keys %stem_ok;
    return @stems[0,1] if @stems >= 2;
    return;
}

sub open_text_reader {
    my ($path) = @_;
    my $win = cygpath_to_win($path);
    if ($win =~ /\.gz$/i || $win =~ /\.bgz$/i || $win =~ /\.bgzip$/i) {
        my $src = $win;
        my $fh;
        my $gzip = external_gzip_path();
        if ($gzip && open($fh, '-|', $gzip, '-dc', $src)) {
            return $fh;
        }
        $fh = IO::Uncompress::Gunzip->new($win)
          or die "Cannot open gzip input $win: $GunzipError\n";
        return $fh;
    }
    open my $fh, '<', $win or die "Cannot read $win: $!\n";
    return $fh;
}

sub external_gzip_path {
    return $EXTERNAL_GZIP if length $EXTERNAL_GZIP;
    my @candidates = grep { length } map {
        my $x = $_;
        $x =~ s/[\r\n]+$//;
        $x;
    } qx(which gzip 2>/dev/null);
    for my $cand (@candidates) {
        next unless -x $cand;
        $EXTERNAL_GZIP = $cand;
        last;
    }
    return $EXTERNAL_GZIP;
}

sub resolve_header_aliases_local {
    my ($cols, $idx, $source_file, $header, $alias_overrides) = @_;
    return resolve_raw_header_aliases(
        cols            => $cols,
        idx             => $idx,
        source_file     => $source_file,
        header          => $header,
        alias_overrides => $alias_overrides,
    );
}

sub normalize_header_local {
    return normalize_header_name(@_);
}

sub find_header_actual {
    my ($cols, $aliases) = @_;
    my %norm_to_actual;
    for my $col (@{$cols}) {
        $norm_to_actual{ normalize_header_local($col) } ||= $col;
    }
    for my $alias (@{$aliases || []}) {
        my $norm = normalize_header_local($alias);
        return $norm_to_actual{$norm} if exists $norm_to_actual{$norm};
    }
    return undef;
}

sub infer_precomputed_spec_from_dir {
    my (%args) = @_;
    my @files = @{ $args{files} || [] };
    my @candidates;
    my $total = scalar @files;
    my $idx = 0;
    for my $path (@files) {
        $idx++;
        progress_note("Inspecting precomputed GWAS [$idx/$total]: " . basename($path));
        my $info = eval { inspect_precomputed_diff_file($path) };
        if ($@) {
            chomp(my $err = $@);
            progress_note("Skipped precomputed candidate " . basename($path) . ": $err");
            next;
        }
        push @candidates, { path => $path, %{$info} };
        progress_note(
            "Accepted precomputed candidate " . basename($path)
            . " as " . $info->{mode}
        );
    }
    return undef unless @candidates;

    @candidates = sort {
        ($b->{mode} eq 'precomputed_diff_stdized' ? 1 : 0) <=> ($a->{mode} eq 'precomputed_diff_stdized' ? 1 : 0)
        || length($b->{path}) <=> length($a->{path})
    } @candidates;
    my $best = $candidates[0];

    my @pairs = build_pairs_from_pair_tags($best->{sample_pairs}, $args{study_prefix});
    die "Could not infer pair definitions from precomputed file $best->{path}\n" unless @pairs;
    my @prefixes = map { $_->{prefix} } @pairs;

    my $project_tag = safe_name($args{project_tag_override} || infer_project_tag_from_dir($args{gwas_dir}, []));
    my $artifact_stem = safe_name($args{artifact_stem_override} || ($project_tag . '_effects'));
    my $detected_build = detect_reference_build_profile(
        file_paths    => [ map { $_->{path} } @candidates ],
        header_lines  => [ map { $_->{header} } @candidates ],
        default_build => 'hg38',
    );
    my $spec = {
        source_mode           => $best->{mode},
        project_tag           => $project_tag,
        artifact_stem         => $artifact_stem,
        reference_build       => ($detected_build->{build} || 'hg38'),
        input_dir             => $args{gwas_dir},
        output_dir            => $args{gwas_dir},
        workdir               => $args{workdir},
        cygwin_bash           => '/bin/bash',
        threshold             => 0.05,
        keep_remote_plot_data => 1,
        top_hit_focus_prefix  => $prefixes[0],
        top_hit_signal_thrshd => '1e-6',
        top_hit_dist_bp       => '1e6',
        local_window_bp       => '1e7',
        local_gtf_window_bp   => '1e7',
        include_non_protein_coding_genes_in_local_gtf => 0,
        open_result           => 0,
        clean_oda_input       => 1,
        pair_col              => 'PAIR_TAG',
        base_cols             => [ qw(CHR BP A1 A2 SNP) ],
        value_fields          => [
            map(($_ . '_BETA'), @{ $best->{side_stems} }),
            'DIFF_BETA',
            map(($_ . '_SE'), @{ $best->{side_stems} }),
            'DIFF_SE',
            map(($_ . '_P'), @{ $best->{side_stems} }),
            'DIFF_P',
            ($best->{mode} eq 'precomputed_diff_stdized' ? ('STD_DIFF_Z', 'STD_DIFF_P') : ()),
        ],
        filter_fields         => [
            map(($_ . '_P'), @{ $best->{side_stems} }),
            'DIFF_P',
            ($best->{mode} eq 'precomputed_diff_stdized' ? ('STD_DIFF_P') : ()),
        ],
        pairs                 => \@pairs,
        ($best->{mode} eq 'precomputed_diff_stdized'
            ? (input_stdized => $best->{path})
            : (input_diff => $best->{path})),
        auto_detect_summary   => {
            generated_from_dir => $args{gwas_dir},
            selected_file      => basename($best->{path}),
            selected_mode      => $best->{mode},
            sample_pair_tags   => $best->{sample_pairs},
            side_stems         => $best->{side_stems},
        },
    };
    my $spec_path = normalize_unix_path($args{spec_out} || "$args{configs_dir}/auto_${artifact_stem}_from_dir.spec.json");
    progress_note(
        "Selected precomputed diff source: " . basename($best->{path})
        . " (" . $best->{mode} . "); pairs=" . scalar(@pairs)
        . ", output=$spec_path"
    );
    return ($spec, $spec_path);
}

sub progress_note {
    my ($msg) = @_;
    my $stamp = strftime('%Y-%m-%d %H:%M:%S', localtime);
    print STDERR "[auto_prepare_and_run_diff_gwas $stamp] $msg\n";
}

sub verify_diff_gwas_deps {
    my ($deps_dir) = @_;
    my @required = (
        'DiffGWASConfig.pm',
        'merge_pgc_vcf_sumstats_long.pl',
        'convert_merged_gwas_to_plotwide.pl',
        'diff_pairwise_gwas.pl',
        'standardize_diff_gwas_zscore.pl',
        'extract_significant_diff_gwas.pl',
        'extract_single_snp_wide_diff_gwas.pl',
        'emit_diff_gwas_runner_env.pl',
        'generate_sas_wide_import_include.pl',
        'render_sas_template.pl',
        'sort_long_gwas_by_coord.sh',
        'run_sas_oda_manhattan4diffgwas_download_png.sh',
        'run_sas_oda_manhattan4diffgwas_uncompressed_download_png.sh',
        'run_sas_oda_local_top_hits_manhattan_download_png.sh',
        'run_sas_oda_local_top_hits_with_gtf_download_html.sh',
        'run_sas_oda_single_snp_with_gtf_download_html.sh',
        'SAS_ODA_Runner.pm',
        'run_sas_oda_manhattan4diffgwas.sas',
        'run_sas_oda_manhattan4diffgwas_uncompressed.sas',
        'run_sas_oda_local_top_hits_manhattan.sas',
        'run_sas_oda_local_top_hits_with_gtf.sas',
        'run_sas_oda_single_snp_with_gtf.sas',
        'Manhattan4DiffGWASs_png.sas',
        'Lattice_gscatter_over_bed_track.sas',
    );
    die "Required dependency directory not found: $deps_dir\n"
      unless -d cygpath_to_win($deps_dir);
    my @missing;
    for my $rel (@required) {
        next unless defined $rel && length $rel;
        my $path = normalize_unix_path(File::Spec->catfile($deps_dir, $rel));
        push @missing, $rel unless -e cygpath_to_win($path);
    }
    my $oda_helper = resolve_oda_helper_unix(script_root_dir());
    push @missing, 'run_sas_codes_or_script_in_ODA.pl'
      unless defined $oda_helper && length $oda_helper && -e cygpath_to_win($oda_helper);
    return if @missing == 0;
    die "Missing required DiffGWASDeps files under $deps_dir:\n  - "
      . join("\n  - ", @missing)
      . "\nPlease restore the missing files before running auto_prepare_and_run_diff_gwas.pl\n";
}

sub build_pairs_from_pair_tags {
    my ($pair_tags, $study_prefix) = @_;
    my @pairs;
    my %prefix_seen;
    for my $pair_tag (@{ $pair_tags || [] }) {
        my ($g1, $g2) = parse_pair_tag($pair_tag);
        next unless length $g1 && length $g2;
        my $code1 = infer_group_code_from_text($g1);
        my $code2 = infer_group_code_from_text($g2);
        my $tag1 = compose_group_tag($study_prefix, $code1);
        my $tag2 = compose_group_tag($study_prefix, $code2);
        my $label1 = expand_group_code($code1);
        my $label2 = expand_group_code($code2);
        my $prefix = uc(safe_name($label1) . '_' . safe_name($label2));
        my $base = $prefix;
        my $suffix = 2;
        while ($prefix_seen{$prefix}) {
            $prefix = $base . '_' . $suffix++;
        }
        $prefix_seen{$prefix} = 1;
        push @pairs, {
            pair_tag     => $pair_tag,
            group1       => $tag1,
            group2       => $tag2,
            group1_label => $label1,
            group2_label => $label2,
            prefix       => $prefix,
            label        => $label1 . ' vs ' . $label2,
        };
    }
    return @pairs;
}

sub parse_pair_tag {
    my ($tag) = @_;
    return ('','') unless defined $tag && length $tag;
    if ($tag =~ /(.+?)_vs_(.+)/i) {
        return ($1, $2);
    }
    my @parts = split /[|:\/-]+/, $tag;
    return @parts >= 2 ? ($parts[0], $parts[1]) : ('','');
}

sub infer_group_tag_from_file {
    my ($path, $study_prefix) = @_;
    my $base = lc basename($path);
    my $code = infer_group_code_from_text($base);
    return compose_group_tag($study_prefix, $code);
}

sub infer_group_label_from_file {
    my ($path, $tag) = @_;
    my $base = lc basename($path);
    my $code = infer_group_code_from_text($base);
    return expand_group_code($code) || ($tag =~ s/^GWAS_//r);
}

sub infer_group_code_from_text {
    my ($text) = @_;
    my @rules = (
        [ qr/female/, 'FEMALE' ],
        [ qr/\bfem\b|_fem_|female/, 'FEMALE' ],
        [ qr/\bmal\b|_mal_|male/, 'MALE' ],
        [ qr/afram|african|afr\b/, 'AFR' ],
        [ qr/asian|eastasian|asn\b/, 'ASN' ],
        [ qr/european|eur\b/, 'EUR' ],
        [ qr/latino|latam|lat\b|amr\b/, 'LAT' ],
        [ qr/primary/, 'PRIMARY' ],
        [ qr/core/, 'CORE' ],
        [ qr/allsex|all_sex|bothsex|combined/, 'ALL' ],
        [ qr/ukbb/, 'UKBB' ],
    );
    for my $rule (@rules) {
        return $rule->[1] if $text =~ $rule->[0];
    }
    my $stem = $text;
    $stem =~ s/\.(?:tsv|txt|sumstats|assoc|vcf)(?:\.(?:gz|bgz|bgzip))?$//i;
    $stem =~ s/\.vcf$//i;
    my @tokens = grep { length && $_ !~ /^(?:pgc\d*|scz|wave\d+|public|autosome|chr\d+|v\d+|tsv|txt|gz|vcf)$/i }
      split /[^A-Za-z0-9]+/, $stem;
    my $token = @tokens ? uc($tokens[0]) : 'GROUP';
    $token =~ s/[^A-Z0-9]+/_/g;
    return $token;
}

sub expand_group_code {
    my ($code) = @_;
    my %map = (
        AFR     => 'AFR',
        ASN     => 'ASN',
        EUR     => 'EUR',
        LAT     => 'LAT',
        FEMALE  => 'Female',
        MALE    => 'Male',
        ALL     => 'All',
        PRIMARY => 'Primary',
        CORE    => 'Core',
    );
    return $map{$code} || $code;
}

sub infer_project_tag_from_dir {
    my ($dir, $groups) = @_;
    my $base = basename($dir);
    return safe_name($base . '_diff');
}

sub select_compatible_groups {
    my ($groups, $warnings) = @_;
    my @groups = @{ $groups || [] };
    return @groups unless @groups >= 2;

    my %by_signature;
    for my $g (@groups) {
        my $sig = join('|', ($g->{format_class} // 'UNKNOWN'), ($g->{effect_metric} // 'UNKNOWN'));
        push @{ $by_signature{$sig} }, $g;
    }
    my @sigs = sort {
        scalar(@{ $by_signature{$b} }) <=> scalar(@{ $by_signature{$a} })
          || $a cmp $b
    } keys %by_signature;
    return @groups if @sigs <= 1;

    my $best_sig = $sigs[0];
    my @best = @{ $by_signature{$best_sig} };
    if (@best >= 2 && @best < @groups) {
        push @{ $warnings || [] },
          'Auto-selection restricted the raw GWAS set to the compatible subgroup '
          . join(', ', map { $_->{tag} } @best)
          . " because the directory mixed multiple schema/effect families. Selected signature: $best_sig.";
        return @best;
    }
    return @groups;
}

sub infer_study_prefix_from_dir {
    my ($dir) = @_;
    my $base = lc basename($dir);
    return 'SCZ_W3' if $base =~ /scz/ && $base =~ /ancestry|sex/;
    return 'GWAS';
}

sub compose_group_tag {
    my ($study_prefix, $code) = @_;
    $study_prefix ||= 'GWAS';
    $code ||= 'GROUP';
    return $study_prefix . '_' . $code;
}

sub win_to_cygpath {
    my ($path) = @_;
    return '' unless defined $path;
    return normalize_unix_path($path) if $path =~ m{^/mnt/};
    if ($path =~ /^([A-Za-z]):[\\\/]?(.*)$/) {
        my ($drive, $rest) = (lc($1), $2 // '');
        $rest =~ s{\\}{/}g;
        return "/mnt/$drive/$rest";
    }
    return normalize_unix_path($path);
}

sub build_pair_info {
    my ($pairs) = @_;
    die "pairs must be a non-empty array\n" unless ref($pairs) eq 'ARRAY' && @{$pairs};
    my (%prefix_seen, %pair_map, %group_rep);
    my (@prefix_order, @labels, @gtf_labels, @pair_defs);
    for my $pair (@{$pairs}) {
        die "Each pair must be an object\n" unless ref($pair) eq 'HASH';
        my $pair_tag = $pair->{pair_tag} // '';
        my $group1 = $pair->{group1} // '';
        my $group2 = $pair->{group2} // '';
        my $prefix = $pair->{prefix} // '';
        my $label = $pair->{label} // $prefix;
        my $group1_label = $pair->{group1_label};
        my $group2_label = $pair->{group2_label};
        $group1_label = infer_group_label($group1) unless defined $group1_label && length $group1_label;
        $group2_label = infer_group_label($group2) unless defined $group2_label && length $group2_label;
        die "pair_tag/group1/group2/prefix are required for each pair\n"
          unless length $pair_tag && length $group1 && length $group2 && length $prefix;
        die "Duplicate prefix in pairs: $prefix\n" if $prefix_seen{$prefix}++;
        $pair_map{$pair_tag} = $prefix;
        push @prefix_order, $prefix;
        push @labels, $label;
        my $gtf = $pair->{gtf_label};
        $gtf = $label unless defined $gtf && length $gtf;
        $gtf =~ s/\s+/_/g;
        push @gtf_labels, $gtf;
        push @pair_defs, {
            %{$pair},
            group1_label => $group1_label,
            group2_label => $group2_label,
        };
        $group_rep{$group1} ||= {
            pvar  => "${prefix}_GROUP1_P",
            betavar => "${prefix}_GROUP1_BETA",
            zvar  => "${prefix}_GROUP1_Z",
            frq_a_var => "${prefix}_GROUP1_FRQ_A",
            frq_u_var => "${prefix}_GROUP1_FRQ_U",
            label => $group1_label,
        };
        $group_rep{$group2} ||= {
            pvar  => "${prefix}_GROUP2_P",
            betavar => "${prefix}_GROUP2_BETA",
            zvar  => "${prefix}_GROUP2_Z",
            frq_a_var => "${prefix}_GROUP2_FRQ_A",
            frq_u_var => "${prefix}_GROUP2_FRQ_U",
            label => $group2_label,
        };
    }
    return {
        pair_map     => \%pair_map,
        prefix_order => \@prefix_order,
        labels       => \@labels,
        gtf_labels   => \@gtf_labels,
        pair_defs    => \@pair_defs,
        group_rep    => \%group_rep,
    };
}

sub parse_display_gwas_list {
    my ($raw) = @_;
    return () unless defined $raw && length $raw;
    $raw =~ s/[|;]/,/g;
    return grep { length } map {
        my $token = $_ // '';
        $token =~ s/^\s+|\s+$//g;
        $token;
    } split /,/, $raw;
}

sub normalize_display_gwas_token {
    my ($token) = @_;
    $token = '' unless defined $token;
    $token =~ s/^\s+|\s+$//g;
    $token =~ s/[[:space:]]+/_/g;
    $token = safe_name($token);
    return uc($token);
}

sub add_display_track_aliases {
    my ($lookup, $entry, @aliases) = @_;
    for my $alias (@aliases) {
        next unless defined $alias && length $alias;
        my $norm = normalize_display_gwas_token($alias);
        next unless length $norm;
        $lookup->{$norm} ||= $entry;
    }
}

sub build_display_track_catalog {
    my ($pair_info, $spec) = @_;
    my @catalog;
    my %lookup;
    my @prefixes = @{ $pair_info->{prefix_order} || [] };
    my @labels = @{ $pair_info->{labels} || [] };
    my @gtf_labels = @{ $pair_info->{gtf_labels} || [] };

    for my $i (0 .. $#prefixes) {
        my $prefix = $prefixes[$i];
        my $label = $labels[$i] // $prefix;
        my $gtf_label = $gtf_labels[$i] // safe_name($label);
        my $entry = {
            id => $prefix,
            kind => 'std',
            prefix => $prefix,
            pvar => $prefix . '_STD_P',
            zvar => $prefix . '_STD_Z',
            betavar => '',
            manhattan_label => $label . ' standardized diff P',
            gtf_label => $gtf_label,
        };
        push @catalog, $entry;
        add_display_track_aliases(
            \%lookup,
            $entry,
            $entry->{id},
            $prefix . '_DIFF',
            $prefix . '_STD',
            $prefix . '_STD_P',
            $prefix . '_STD_Z',
            $label,
            $label . ' diff',
            $label . ' standardized diff',
            $gtf_label,
        );
    }

    for my $group (sort keys %{ $pair_info->{group_rep} || {} }) {
        my $rep = $pair_info->{group_rep}{$group};
        my $safe_label = safe_name($rep->{label});
        my $entry = {
            id => $safe_label,
            kind => 'group',
            group_key => $group,
            prefix => '',
            pvar => $safe_label . '_P',
            zvar => $safe_label . '_Z',
            betavar => $safe_label . '_BETA',
            manhattan_label => $rep->{label} . ' association P',
            gtf_label => $safe_label,
        };
        push @catalog, $entry;
        add_display_track_aliases(
            \%lookup,
            $entry,
            $entry->{id},
            $group,
            $safe_label . '_P',
            $safe_label . '_Z',
            $rep->{label},
            $rep->{label} . ' association',
        );
    }

    for my $track (@{ cfg_or($spec || {}, 'extra_tracks', []) || [] }) {
        next unless ref($track) eq 'HASH';
        my $id = $track->{id} || next;
        my $label = $track->{label} || $id;
        my $entry = {
            id => $id,
            kind => ($track->{kind} || 'extra'),
            group_key => '',
            prefix => '',
            pvar => $track->{pvar},
            zvar => $track->{zvar},
            betavar => ($track->{betavar} || ''),
            sevar => ($track->{sevar} || ''),
            manhattan_label => ($track->{manhattan_label} || ($label . ' association P')),
            gtf_label => ($track->{gtf_label} || safe_name($label)),
            label => $label,
        };
        push @catalog, $entry;
        add_display_track_aliases(
            \%lookup,
            $entry,
            $entry->{id},
            $entry->{pvar},
            $entry->{zvar},
            $label,
            $entry->{gtf_label},
        );
    }

    return (\@catalog, \%lookup);
}

sub resolve_display_track_selection {
    my (%args) = @_;
    my $spec = $args{spec};
    my $pair_info = $args{pair_info};
    my $override = $args{display_gwas_override} // '';
    my $raw = length($override) ? $override : cfg_or($spec, 'display_gwas', '');
    my ($catalog, $lookup) = build_display_track_catalog($pair_info, $spec);
    my @tokens = parse_display_gwas_list($raw);

    return {
        explicit => 0,
        raw => '',
        tracks => [ @{$catalog} ],
        ids => [ map { $_->{id} } @{$catalog} ],
        canonical => join(',', map { $_->{id} } @{$catalog}),
        available => [ map { $_->{id} } @{$catalog} ],
    } unless @tokens;

    my @selected;
    my %seen;
    for my $token (@tokens) {
        my $norm = normalize_display_gwas_token($token);
        my $entry = $lookup->{$norm};
        die "Unknown --display-gwas value '$token'. Available values: "
          . join(', ', map { $_->{id} } @{$catalog}) . "\n"
          unless $entry;
        next if $seen{ $entry->{id} }++;
        push @selected, { %{$entry} };
    }

    die "No valid display tracks remained after applying --display-gwas\n" unless @selected;

    return {
        explicit => 1,
        raw => $raw,
        tracks => \@selected,
        ids => [ map { $_->{id} } @selected ],
        canonical => join(',', map { $_->{id} } @selected),
        available => [ map { $_->{id} } @{$catalog} ],
    };
}

sub display_selection_variant_tag {
    my ($selection) = @_;
    return '' unless ref($selection) eq 'HASH' && $selection->{explicit};
    my @ids = @{ $selection->{ids} || [] };
    return '' unless @ids;
    my $slug = safe_name(join('_', @ids));
    $slug = substr($slug, 0, 48) if length($slug) > 48;
    my $hash = substr(md5_hex(join(',', @ids)), 0, 10);
    return join('_', grep { length } ('display', $slug, $hash));
}

sub build_forest_group_tracks {
    my (%args) = @_;
    my $selection = $args{selection} || {};
    my $pair_info = $args{pair_info} || {};
    my %group_rep = %{ $pair_info->{group_rep} || {} };
    my @pair_defs = @{ $pair_info->{pair_defs} || [] };
    my @tracks = @{ $selection->{tracks} || [] };

    my %pair_by_prefix = map { normalize_display_gwas_token($_->{prefix}) => $_ } @pair_defs;
    my @group_order;
    my %seen;

    if ($selection->{explicit}) {
        for my $track (@tracks) {
            next unless ref($track) eq 'HASH';
            if (($track->{kind} || '') eq 'group') {
                my $group_key = $track->{group_key} || '';
                next unless length($group_key) && exists $group_rep{$group_key};
                next if $seen{$group_key}++;
                push @group_order, $group_key;
                next;
            }
            if (($track->{kind} || '') eq 'std') {
                my $pair = $pair_by_prefix{ normalize_display_gwas_token($track->{prefix} || $track->{id} || '') };
                next unless $pair;
                for my $group_key ($pair->{group1}, $pair->{group2}) {
                    next unless defined $group_key && length($group_key) && exists $group_rep{$group_key};
                    next if $seen{$group_key}++;
                    push @group_order, $group_key;
                }
            }
        }
    }

    if (!@group_order) {
        for my $pair (@pair_defs) {
            for my $group_key ($pair->{group1}, $pair->{group2}) {
                next unless defined $group_key && length($group_key) && exists $group_rep{$group_key};
                next if $seen{$group_key}++;
                push @group_order, $group_key;
            }
        }
    }

    my @forest_tracks;
    for my $group_key (@group_order) {
        my $rep = $group_rep{$group_key} || next;
        my $betavar = $rep->{betavar} || '';
        my $sevar = $betavar;
        $sevar =~ s/_BETA$/_SE/ if length $sevar;
        push @forest_tracks, {
            id        => safe_name($rep->{label}),
            group_key => $group_key,
            label     => $rep->{label},
            pvar      => $rep->{pvar},
            betavar   => $betavar,
            sevar     => $sevar,
            zvar      => $rep->{zvar},
        };
    }
    return \@forest_tracks;
}

sub append_variant_to_filename {
    my ($name, $variant) = @_;
    return $name unless defined $variant && length $variant;
    return '' unless defined $name;
    if ($name =~ /^(.*?)(\.[^.]+)$/) {
        return $1 . '_' . $variant . $2;
    }
    return $name . '_' . $variant;
}

sub infer_effect_metric_label_from_vars {
    my (@vars) = @_;
    @vars = grep { defined $_ && length $_ } @vars;
    return 'Effect metric' unless @vars;
    my $all_z = 1;
    my $all_beta = 1;
    my $all_or = 1;
    for my $var (@vars) {
        my $u = uc($var);
        $all_z = 0 unless $u =~ /(?:^|_)(?:Z|ZSCORE)(?:$|_)/;
        $all_beta = 0 unless $u =~ /(?:^|_)(?:BETA|EFFECT)(?:$|_)/;
        $all_or = 0 unless $u =~ /(?:^|_)(?:OR|ODDSRATIO)(?:$|_)/;
    }
    return 'Z score' if $all_z;
    return 'Beta' if $all_beta;
    return 'Odds ratio' if $all_or;
    return 'Effect metric';
}

sub build_merge_config {
    my ($spec, $generated) = @_;
    return {
        input_dir => normalize_unix_path($spec->{input_dir}),
        output    => $generated->{merge_output},
        manifest  => $generated->{merge_manifest},
        (ref($spec->{raw_column_aliases}) eq 'HASH' && keys %{ $spec->{raw_column_aliases} || {} }
            ? (raw_column_aliases => $spec->{raw_column_aliases})
            : ()),
        groups    => $spec->{groups},
    };
}

sub build_diff_config {
    my ($spec, $generated) = @_;
    my %pairs = map {
        $_->{pair_tag} => [ $_->{group1}, $_->{group2} ]
    } @{ $spec->{pairs} };
    return {
        input         => $generated->{sorted_output},
        output        => $generated->{diff_output},
        manifest      => $generated->{diff_manifest},
        base_cols     => [ qw(CHR BP A1 A2 SNP) ],
        group_tag_col => 'GWAS_TAG',
        beta_col      => 'BETA',
        se_col        => 'SE',
        p_col         => 'P',
        rho           => cfg_or($spec, 'rho', 0),
        pairs         => \%pairs,
    };
}

sub build_preset_config {
    my ($spec, $generated, $pair_info, $threshold, $window_bp) = @_;
    my $source_mode = cfg_or($spec, 'source_mode', 'raw_pgc_vcf_sumstats');
    if ($source_mode eq 'merged_gwas_table') {
        my @wide_columns = build_merged_wide_columns($spec, $pair_info);
        return {
            project_tag   => cfg_or($spec, 'artifact_stem', cfg_or($spec, 'project_tag', 'diff_gwas')),
            input         => normalize_unix_path(cfg_or($spec, 'input_merged', '')),
            output        => $generated->{wide_output},
            manifest      => $generated->{wide_manifest},
            output_dir    => normalize_unix_path(cfg_or($spec, 'workdir', script_root_dir())),
            threshold     => $threshold + 0,
            window_bp     => $window_bp + 0,
            base_cols     => [ qw(CHR BP A1 A2 SNP) ],
            wide_columns  => \@wide_columns,
            char_lengths  => {
                A1  => 8,
                A2  => 8,
                SNP => 128,
            },
            alias_map      => build_alias_map($pair_info),
            post_alias_map => build_post_alias_map($pair_info),
            pair_map       => $pair_info->{pair_map},
            prefix_order   => $pair_info->{prefix_order},
            merged_base_cols   => cfg_or($spec, 'merged_base_cols', {}),
            merged_group_tracks => cfg_or($spec, 'merged_group_tracks', []),
            merged_extra_tracks => cfg_or($spec, 'merged_extra_tracks', []),
            pairs          => cfg_or($spec, 'pairs', []),
            extra_tracks   => cfg_or($spec, 'extra_tracks', []),
        };
    }
    my ($value_fields, $filter_fields) = detect_stdized_value_and_filter_fields($generated->{stdized_output});
    return {
        project_tag   => cfg_or($spec, 'artifact_stem', cfg_or($spec, 'project_tag', 'diff_gwas')),
        input         => $generated->{stdized_output},
        output        => $generated->{wide_output},
        manifest      => $generated->{wide_manifest},
        output_dir    => normalize_unix_path(cfg_or($spec, 'workdir', '/mnt/g/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper')),
        threshold     => $threshold + 0,
        window_bp     => $window_bp + 0,
        pair_col      => cfg_or($spec, 'pair_col', 'PAIR_TAG'),
        base_cols     => cfg_or($spec, 'base_cols', [ qw(CHR BP A1 A2 SNP) ]),
        value_fields  => cfg_or($spec, 'value_fields', $value_fields),
        filter_fields => cfg_or($spec, 'filter_fields', $filter_fields),
        char_lengths  => {
            A1  => 8,
            A2  => 8,
            SNP => 40,
        },
        alias_map      => build_alias_map($pair_info),
        post_alias_map => build_post_alias_map($pair_info),
        pair_map       => $pair_info->{pair_map},
        prefix_order   => $pair_info->{prefix_order},
    };
}

sub build_merged_wide_columns {
    my ($spec, $pair_info) = @_;
    my @cols = qw(CHR BP A1 A2 SNP);
    for my $prefix (@{ $pair_info->{prefix_order} || [] }) {
        push @cols,
          "${prefix}_GROUP1_BETA",
          "${prefix}_GROUP1_SE",
          "${prefix}_GROUP1_P",
          "${prefix}_GROUP1_Z",
          "${prefix}_GROUP2_BETA",
          "${prefix}_GROUP2_SE",
          "${prefix}_GROUP2_P",
          "${prefix}_GROUP2_Z",
          "${prefix}_DIFF_BETA",
          "${prefix}_DIFF_SE",
          "${prefix}_DIFF_P",
          "${prefix}_STD_DIFF_Z",
          "${prefix}_STD_DIFF_P";
    }
    for my $track (@{ cfg_or($spec, 'extra_tracks', []) || [] }) {
        next unless ref($track) eq 'HASH';
        my $id = $track->{id} || next;
        push @cols, "${id}_BETA" if defined $track->{betavar} && length($track->{betavar});
        push @cols, "${id}_SE"   if defined $track->{sevar}   && length($track->{sevar});
        push @cols, "${id}_P"    if defined $track->{pvar}    && length($track->{pvar});
        push @cols, "${id}_Z"    if defined $track->{zvar}    && length($track->{zvar});
    }
    my %seen;
    return grep { !$seen{$_}++ } @cols;
}

sub detect_stdized_value_and_filter_fields {
    my ($stdized_output) = @_;
    my @default_value_fields = qw(
      GROUP1_BETA GROUP2_BETA DIFF_BETA
      GROUP1_SE GROUP2_SE DIFF_SE
      GROUP1_P GROUP2_P DIFF_P STD_DIFF_Z STD_DIFF_P
      GROUP1_FRQ_A GROUP1_FRQ_U GROUP2_FRQ_A GROUP2_FRQ_U
      GROUP1_INFO GROUP2_INFO
    );
    my @default_filter_fields = qw(GROUP1_P GROUP2_P DIFF_P STD_DIFF_P);
    return (\@default_value_fields, \@default_filter_fields)
      unless defined $stdized_output && length $stdized_output;

    my %idx;
    eval {
        my $fh = open_text_reader($stdized_output);
        my $header = <$fh>;
        die "Empty standardized diff file: $stdized_output\n" unless defined $header;
        chomp $header;
        $header =~ s/\r$//;
        my @cols = split /\t/, $header, -1;
        %idx = map { $cols[$_] => 1 } 0 .. $#cols;
        close $fh;
        1;
    } or do {
        warn "Unable to inspect standardized header for $stdized_output: $@";
        return (\@default_value_fields, \@default_filter_fields);
    };

    my @families = (
        {
            required_value_fields => [ qw(GROUP1_BETA GROUP2_BETA DIFF_BETA GROUP1_SE GROUP2_SE DIFF_SE GROUP1_P GROUP2_P DIFF_P STD_DIFF_Z STD_DIFF_P) ],
            optional_value_fields => [ qw(GROUP1_FRQ_A GROUP1_FRQ_U GROUP2_FRQ_A GROUP2_FRQ_U GROUP1_INFO GROUP2_INFO) ],
            filter_fields => [ qw(GROUP1_P GROUP2_P DIFF_P STD_DIFF_P) ],
        },
        {
            required_value_fields => [ qw(FEMALE_BETA MALE_BETA DIFF_BETA FEMALE_SE MALE_SE DIFF_SE FEMALE_P MALE_P DIFF_P STD_DIFF_Z STD_DIFF_P) ],
            optional_value_fields => [ qw(FEMALE_FRQ_A FEMALE_FRQ_U MALE_FRQ_A MALE_FRQ_U FEMALE_INFO MALE_INFO) ],
            filter_fields => [ qw(FEMALE_P MALE_P DIFF_P STD_DIFF_P) ],
        },
    );
    for my $family (@families) {
        my $ok = 1;
        for my $col (@{ $family->{required_value_fields} }) {
            if (!exists $idx{$col}) {
                $ok = 0;
                last;
            }
        }
        if ($ok) {
            my @value_fields = (
                @{ $family->{required_value_fields} },
                grep { exists $idx{$_} } @{ $family->{optional_value_fields} || [] },
            );
            return (\@value_fields, $family->{filter_fields});
        }
    }
    return (\@default_value_fields, \@default_filter_fields);
}

sub build_runner_config {
    my (%args) = @_;
    my $spec = $args{spec};
    my $source_mode = cfg_or($spec, 'source_mode', 'raw_pgc_vcf_sumstats');
    my $generated = $args{generated};
    my $pair_info = $args{pair_info};
    my $project_tag = $args{project_tag};
    my $reference_build_profile = $args{reference_build_profile} || {};
    my $focus_prefix = $args{focus_prefix};
    my $top_hit_threshold = $args{top_hit_threshold};
    my $top_hit_dist_bp = $args{top_hit_dist_bp};
    my $top_hit_max_loci = $args{top_hit_max_loci};
    my $local_window_bp = $args{local_window_bp};
    my $local_gtf_window_bp = $args{local_gtf_window_bp};
    my $local_manhattan_angle4xaxis_label_override = $args{local_manhattan_angle4xaxis_label_override};
    my $local_manhattan_xgrp_y_pos_override = $args{local_manhattan_xgrp_y_pos_override};
    my $local_manhattan_yoffset_top_override = $args{local_manhattan_yoffset_top_override};
    my $local_manhattan_yoffset_bottom_override = $args{local_manhattan_yoffset_bottom_override};
    my $local_manhattan_fontsize_override = $args{local_manhattan_fontsize_override};
    my $local_manhattan_y_axis_label_size_override = $args{local_manhattan_y_axis_label_size_override};
    my $local_manhattan_y_axis_value_size_override = $args{local_manhattan_y_axis_value_size_override};
    my $target_snps_override = $args{target_snps_override} // '';
    my $target_snp_genes_override = $args{target_snp_genes_override} // '';
    my $display_gwas_override = $args{display_gwas_override} // '';
    my $local_gtf_label_snps_override = $args{local_gtf_label_snps_override} // '';
    my $local_gtf_label_layout_override = $args{local_gtf_label_layout_override} // '';
    my $local_gtf_yaxis_offset4max_override = $args{local_gtf_yaxis_offset4max_override} // '';
    my $local_gtf_yoffset4textlabels_override = $args{local_gtf_yoffset4textlabels_override} // '';
    my $get_common_associations = $args{get_common_associations} ? 1 : 0;
    my $common_assoc_top_hit_threshold_override = $args{common_assoc_top_hit_threshold_override} // '';
    my $selection = resolve_display_track_selection(
        spec => $spec,
        pair_info => $pair_info,
        display_gwas_override => $display_gwas_override,
    );
    my @selected_tracks = @{ $selection->{tracks} || [] };
    my @selected_std_tracks = grep { $_->{kind} eq 'std' } @selected_tracks;
    my @selected_group_tracks = grep { $_->{kind} eq 'group' } @selected_tracks;
    my @forest_tracks = @{ build_forest_group_tracks(
        selection => $selection,
        pair_info => $pair_info,
    ) || [] };

    my @prefixes = map { $_->{prefix} } @selected_std_tracks;
    my @labels = map { $_->{manhattan_label} } @selected_tracks;
    my @gtf_labels = map { $_->{gtf_label} } @selected_tracks;
    my @std_pvars = map { $_->{pvar} } @selected_std_tracks;
    my @std_zvars = map { $_->{zvar} } @selected_std_tracks;
    my @group_pvars = map { $_->{pvar} } @selected_group_tracks;
    my @group_betavars = map { $_->{betavar} } @selected_group_tracks;
    my @group_zvars = map { $_->{zvar} } @selected_group_tracks;
    my $track_count = scalar(@selected_tracks) || 1;
    my $local_max_hits_per_fig = $args{local_max_hits_per_fig_override}
      ? $args{local_max_hits_per_fig_override}
      : cfg_or($spec, 'local_max_hits_per_fig', 15);
    $local_max_hits_per_fig = 15 unless defined $local_max_hits_per_fig && $local_max_hits_per_fig =~ /^\d+$/ && $local_max_hits_per_fig > 0;
    $local_max_hits_per_fig = 15 if $local_max_hits_per_fig > 15;
    # Local GTF plotting is much more sensitive to SAS ODA WORK-space pressure
    # than local Manhattan. Enforce one locus per batch/figure so each top hit
    # runs in isolation.
    my $local_gtf_max_hits_per_fig = 1;
    my $manhattan_fig_height = cfg_or($spec, 'manhattan_fig_height', suggest_genomewide_fig_height($track_count));
    my $manhattan_fig_width = cfg_or($spec, 'manhattan_fig_width', 1800);
    my $local_manhattan_fig_height = cfg_or($spec, 'local_manhattan_fig_height', suggest_local_manhattan_fig_height($track_count, $local_max_hits_per_fig));
    my $local_manhattan_fig_width = cfg_or($spec, 'local_manhattan_fig_width', 1800);
    my $common_assoc_nominal_thrshd = cfg_or($spec, 'common_assoc_nominal_thrshd', '0.05');
    my $common_assoc_thrshds = cfg_or($spec, 'common_assoc_signal_thrshds', [ '5e-8', '1e-6', '1e-5' ]);
    $common_assoc_thrshds = [ split(/[,\s]+/, $common_assoc_thrshds) ] unless ref($common_assoc_thrshds) eq 'ARRAY';
    my @common_assoc_thrshds = grep { defined $_ && length $_ } @{ $common_assoc_thrshds };
    @common_assoc_thrshds = ('5e-8', '1e-6', '1e-5') unless @common_assoc_thrshds;
    if (defined $common_assoc_top_hit_threshold_override && length $common_assoc_top_hit_threshold_override) {
        my %seen;
        @common_assoc_thrshds = grep { !$seen{$_}++ } ($common_assoc_top_hit_threshold_override, @common_assoc_thrshds);
    }

    my @mh_labels = map { $_->{manhattan_label} } @selected_tracks;

    my @focus_tracks = @selected_std_tracks ? @selected_std_tracks : @selected_group_tracks;
    die "No selectable Manhattan/GTF tracks were resolved for plotting.\n" unless @focus_tracks;

    my $primary_focus_track;
    if ($selection->{explicit}) {
        $primary_focus_track = $focus_tracks[0];
    }
    elsif (@selected_std_tracks) {
        ($primary_focus_track) = grep { $_->{prefix} eq $focus_prefix } @selected_std_tracks;
        $primary_focus_track ||= $selected_std_tracks[0];
    }
    else {
        $primary_focus_track = $selected_group_tracks[0];
    }

    my $top_hit_mode = $get_common_associations
        ? 'common_association'
        : ($source_mode eq 'merged_gwas_table' && @selected_std_tracks && @selected_group_tracks >= 2
            ? 'common_and_differential'
            : (@selected_std_tracks ? 'differential' : 'single_gwas'));
    my $top_hit_focus_pvar = $get_common_associations ? 'COMMON_ASSOC_P' : $primary_focus_track->{pvar};
    my $top_hit_filter_expr;
    my $top_hit_signal_thrshds;
    if ($get_common_associations) {
                die "Common-association mode requires at least two displayed single-GWAS tracks. "
                  . "Use --display-gwas with at least two GWAS association tracks, or omit it.\n"
                  if @selected_group_tracks < 2;
                my @common_repl_exprs;
                for my $driver_track (@selected_group_tracks) {
                        my $driver_p = $driver_track->{pvar};
                        my $driver_z = $driver_track->{zvar};
                        for my $partner_track (@selected_group_tracks) {
                                next if $partner_track->{id} eq $driver_track->{id};
                                my $partner_p = $partner_track->{pvar};
                                my $partner_z = $partner_track->{zvar};
                                my $same_dir = "((($driver_z>0) and ($partner_z>0)) or (($driver_z<0) and ($partner_z<0)))";
                                push @common_repl_exprs,
                                    "((($driver_p>0) and ($driver_p=$top_hit_focus_pvar)) and (($partner_p>0) and ($partner_p<=$common_assoc_nominal_thrshd)) and $same_dir)";
                        }
                }
                $top_hit_filter_expr = join(' or ', @common_repl_exprs);
                $top_hit_signal_thrshds = join(' ', @common_assoc_thrshds);

                # Historically the runner filter expression sometimes carried a
                # redundant standardized-differential non-significance clause such as
                # "(PREFIX_STD_P>=0.5)" which can suppress valid common-association
                # candidates. Remove any accidental "_STD_P>=0.5" clauses when in
                # common-association mode unless explicitly requested in the spec.
                $top_hit_filter_expr =~ s/\s+and\s+\([A-Za-z0-9_]+_STD_P\s*>=\s*0(?:\.0+)?\)//gi;
                $top_hit_filter_expr =~ s/\s+and\s+\([A-Za-z0-9_]+_STD_P\s*>\s*=\s*0(?:\.0+)?\)//gi; # tolerance for weird spacing
                $top_hit_filter_expr =~ s/\s+and\s+\([A-Za-z0-9_]+_STD_P\s*>\s*0(?:\.0+)?\)//gi; # defensive

                # Clean up accidental duplicated whitespace or leading/trailing connectors
                $top_hit_filter_expr =~ s/^\s+|\s+$//g;
                $top_hit_filter_expr =~ s/\s+or\s+\s+or\s+/ or /gi;
    } else {
        $top_hit_filter_expr = join(' or ', map {
            "(($_>0) and ($_<$top_hit_threshold))"
        } map { $_->{pvar} } @focus_tracks);
        if ($top_hit_mode eq 'differential') {
            $top_hit_signal_thrshds = build_differential_threshold_ladder(
                $top_hit_threshold,
                $top_hit_threshold_fallback,
            );
        } else {
            $top_hit_signal_thrshds = $top_hit_threshold;
        }
    }

    my $display_variant = display_selection_variant_tag($selection);
    my $output_prefix = cfg_or($spec, 'output_prefix', $project_tag . '_SAS_manhattan');
    my $local_output_prefix = cfg_or($spec, 'local_output_prefix', $project_tag . '_SAS_local_top_hits_manhattan');
    my $local_top_hits_csv_basename = cfg_or($spec, 'local_top_hits_csv_basename', $project_tag . '_SAS_local_top_hits_manhattan_top_hits.csv');
    my $output_html_basename = cfg_or($spec, 'output_html_basename', $project_tag . '_SAS_local_top_hits_with_gtf.html');
    my $forest_output_prefix = cfg_or($spec, 'forest_output_prefix', $project_tag . '_SAS_top_hits_forest');
    my $forest_top_hits_csv_basename = cfg_or($spec, 'forest_top_hits_csv_basename', $forest_output_prefix . '_top_hits.csv');
    my $forest_output_html_basename = cfg_or($spec, 'forest_output_html_basename', $forest_output_prefix . '.html');
    my $forest_output_manifest_basename = cfg_or($spec, 'forest_output_manifest_basename', $forest_output_prefix . '.manifest.tsv');
    if (length $display_variant) {
        $output_prefix = append_variant_to_filename($output_prefix, $display_variant);
        $local_output_prefix = append_variant_to_filename($local_output_prefix, $display_variant);
        $local_top_hits_csv_basename = append_variant_to_filename($local_top_hits_csv_basename, $display_variant);
        $output_html_basename = append_variant_to_filename($output_html_basename, $display_variant);
        $forest_output_prefix = append_variant_to_filename($forest_output_prefix, $display_variant);
        $forest_top_hits_csv_basename = append_variant_to_filename($forest_top_hits_csv_basename, $display_variant);
        $forest_output_html_basename = append_variant_to_filename($forest_output_html_basename, $display_variant);
        $forest_output_manifest_basename = append_variant_to_filename($forest_output_manifest_basename, $display_variant);
    }

    my $forest_default_hit_class = $top_hit_mode eq 'common_association'
        ? 'COMMON'
        : ($top_hit_mode eq 'differential' || $top_hit_mode eq 'common_and_differential'
            ? 'DIFFERENTIAL'
            : ($top_hit_mode eq 'single_gwas' ? 'SINGLE_GWAS' : 'CUSTOM_TARGET'));

    return {
        PROJECT_TAG => $project_tag,
        REFERENCE_BUILD => ($reference_build_profile->{build} || 'hg38'),
        REFERENCE_BUILD_SOURCE => ($reference_build_profile->{source} || 'fallback_default'),
        REFERENCE_BUILD_EVIDENCE => ($reference_build_profile->{evidence} || ''),
        DATA_GZ => $generated->{wide_output},
        SOURCE_LONG_GZ => ($source_mode eq 'merged_gwas_table' ? '' : $generated->{stdized_output}),
        EXTRACTOR_CONFIG_JSON => $generated->{preset_config},
        DISPLAY_GWAS => $selection->{canonical},
        DISPLAY_GWAS_MODE => (scalar(@selected_tracks) == 1 ? 'single' : 'multi'),
        DISPLAY_GWAS_AVAILABLE => join('|', @{ $selection->{available} || [] }),
        MANHATTAN_GWAS_MODE => (scalar(@selected_tracks) == 1 ? 'single' : 'multi'),
        MANHATTAN_P_VAR => $selected_tracks[0]{pvar},
        MANHATTAN_OTHER_P_VARS => [
            map { $_->{pvar} } @selected_tracks[1 .. $#selected_tracks],
        ],
        MANHATTAN_FIG_HEIGHT => $manhattan_fig_height,
        MANHATTAN_FIG_WIDTH => $manhattan_fig_width,
        MANHATTAN_GWAS_LABEL_NAMES => join('|', @mh_labels),
        TOP_HIT_MODE => $top_hit_mode,
        TOP_HIT_FOCUS_PVAR => $top_hit_focus_pvar,
        TOP_HIT_FILTER_EXPR => $top_hit_filter_expr,
        TOP_HIT_SIGNAL_THRSHD => $top_hit_threshold,
        TOP_HIT_SIGNAL_THRSHDS => $top_hit_signal_thrshds,
        TOP_HIT_SIGNAL_THRSHD_FALLBACK => (
            $top_hit_mode eq 'differential'
              ? differential_threshold_fallback_value($top_hit_threshold, $top_hit_threshold_fallback)
              : ''
        ),
        TOP_HIT_DIST_BP => $top_hit_dist_bp,
        TOP_HIT_MAX_LOCI => (
            defined($top_hit_max_loci) && $top_hit_max_loci =~ /^\d+$/
              ? $top_hit_max_loci
              : 0
        ),
        TOP_HIT_MAF_THRESHOLD => cfg_or($spec, 'top_hit_maf_threshold', ($source_mode eq 'merged_gwas_table' ? 0 : 0.01)),
        TOP_HIT_GNOMAD_FREQ_FILE => cfg_or($spec, 'gnomad_freq_file', ''),
        TOP_HIT_GNOMAD_POP_MAP => population_map_to_string(cfg_or($spec, 'gnomad_population_map', '')),
        TARGET_SNP_LIST => (
            length($target_snps_override)
              ? $target_snps_override
              : cfg_or($spec, 'target_snps', '')
        ),
        TARGET_SNP_GENES => normalize_target_snp_gene_overrides(
            length($target_snp_genes_override)
              ? $target_snp_genes_override
              : cfg_or($spec, 'target_snp_genes', '')
        ),
        COMMON_ASSOC_P_VARS => join(' ', @group_pvars),
        COMMON_ASSOC_BETA_VARS => join(' ', @group_betavars),
        LOCAL_WINDOW_BP => $local_window_bp,
        LOCAL_GTF_WINDOW_BP => $local_gtf_window_bp,
        GTF_DSD => cfg_or($spec, 'gtf_dsd', ($reference_build_profile->{shared_dsd} || 'FM.GTF_HG38')),
        GTF_LOCAL_DSD => cfg_or($spec, 'gtf_local_dsd', ($reference_build_profile->{local_dsd} || 'gtf_hg38')),
        GTF_GZ_URL => cfg_or($spec, 'gtf_gz_url', ($reference_build_profile->{gtf_url} || 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz')),
        LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES => (
            cfg_or($spec, 'include_non_protein_coding_genes_in_local_gtf', 0) ? 1 : 0
        ),
        LOCAL_MANHATTAN_FIG_HEIGHT => $local_manhattan_fig_height,
        LOCAL_MANHATTAN_FIG_WIDTH => $local_manhattan_fig_width,
        LOCAL_MAX_HITS_PER_FIG => $local_max_hits_per_fig,
        LOCAL_GTF_MAX_HITS_PER_FIG => $local_gtf_max_hits_per_fig,
        LOCAL_MANHATTAN_ANGLE4XAXIS_LABEL => (
            length($local_manhattan_angle4xaxis_label_override // '')
              ? $local_manhattan_angle4xaxis_label_override
              : cfg_or($spec, 'local_manhattan_angle4xaxis_label', '')
        ),
        LOCAL_MANHATTAN_XGRP_Y_POS => (
            length($local_manhattan_xgrp_y_pos_override // '')
              ? $local_manhattan_xgrp_y_pos_override
              : cfg_or($spec, 'local_manhattan_xgrp_y_pos', '')
        ),
        LOCAL_MANHATTAN_YOFFSET_TOP => (
            length($local_manhattan_yoffset_top_override // '')
              ? $local_manhattan_yoffset_top_override
              : cfg_or($spec, 'local_manhattan_yoffset_top', '')
        ),
        LOCAL_MANHATTAN_YOFFSET_BOTTOM => (
            length($local_manhattan_yoffset_bottom_override // '')
              ? $local_manhattan_yoffset_bottom_override
              : cfg_or($spec, 'local_manhattan_yoffset_bottom', '')
        ),
        LOCAL_MANHATTAN_FONTSIZE => (
            length($local_manhattan_fontsize_override // '')
              ? $local_manhattan_fontsize_override
              : cfg_or($spec, 'local_manhattan_fontsize', '')
        ),
        LOCAL_MANHATTAN_Y_AXIS_LABEL_SIZE => (
            length($local_manhattan_y_axis_label_size_override // '')
              ? $local_manhattan_y_axis_label_size_override
              : cfg_or($spec, 'local_manhattan_y_axis_label_size', '')
        ),
        LOCAL_MANHATTAN_Y_AXIS_VALUE_SIZE => (
            length($local_manhattan_y_axis_value_size_override // '')
              ? $local_manhattan_y_axis_value_size_override
              : cfg_or($spec, 'local_manhattan_y_axis_value_size', '')
        ),
        GTF_ASSOC_PVARS => join(' ', map { $_->{pvar} } @selected_tracks),
        GTF_ZSCORE_VARS => join(' ', map { $_->{zvar} } @selected_tracks),
        GTF_LABELS => join(' ', @gtf_labels),
        GTF_DIST2SNP => cfg_or($spec, 'gtf_dist2snp', 100000),
        GTF_YAXIS_LABEL => cfg_or($spec, 'gtf_yaxis_label', '-log10(P)'),
        GTF_COLORBAR_LABEL => cfg_or(
            $spec,
            'gtf_colorbar_label',
            infer_effect_metric_label_from_vars(map { $_->{zvar} } @selected_tracks)
        ),
        GTF_YAXIS_OFFSET4MAX => (
            length($local_gtf_yaxis_offset4max_override)
              ? $local_gtf_yaxis_offset4max_override
              : cfg_or($spec, 'gtf_yaxis_offset4max', '')
        ),
        GTF_YOFFSET4TEXTLABELS => (
            length($local_gtf_yoffset4textlabels_override)
              ? $local_gtf_yoffset4textlabels_override
              : cfg_or($spec, 'gtf_yoffset4textlabels', '')
        ),
        GTF_YOFFSET4MAX_DRAWMARKERSONTOP => cfg_or($spec, 'gtf_yoffset4max_drawmarkersontop', ''),
        GTF_LABEL_SNPS => (
            length($local_gtf_label_snps_override)
              ? $local_gtf_label_snps_override
              : cfg_or($spec, 'local_gtf_label_snps', '')
        ),
        GTF_LABEL_LAYOUT => (
            length($local_gtf_label_layout_override)
              ? $local_gtf_label_layout_override
              : cfg_or($spec, 'local_gtf_label_layout', 'auto')
        ),
        FOREST_TRACK_IDS => join('|', map { $_->{id} } @forest_tracks),
        FOREST_TRACK_LABELS => join('|', map { $_->{label} } @forest_tracks),
        FOREST_TRACK_BETA_VARS => join('|', map { $_->{betavar} } @forest_tracks),
        FOREST_TRACK_SE_VARS => join('|', map { $_->{sevar} } @forest_tracks),
        FOREST_TRACK_P_VARS => join('|', map { $_->{pvar} } @forest_tracks),
        FOREST_TRACK_Z_VARS => join('|', map { $_->{zvar} } @forest_tracks),
        FOREST_TRACK_COUNT => scalar(@forest_tracks),
        FOREST_DEFAULT_HIT_CLASS => $forest_default_hit_class,
        FOREST_FIG_WIDTH => cfg_or($spec, 'forest_fig_width', 900),
        FOREST_FIG_HEIGHT => cfg_or($spec, 'forest_fig_height', ''),
        FOREST_DOTSIZE => cfg_or($spec, 'forest_dotsize', 8),
        FOREST_Y_FONT_SIZE => cfg_or($spec, 'forest_y_font_size', 12),
        FOREST_MIN_AXIS => cfg_or($spec, 'forest_min_axis', 0.4),
        FOREST_MAX_AXIS => cfg_or($spec, 'forest_max_axis', 1.6),
        FOREST_XAXIS_VALUE_RANGE => cfg_or($spec, 'forest_xaxis_value_range', '0.4 to 1.6 by 0.2'),
        OUTPUT_PREFIX => $output_prefix,
        HTML_TITLE => cfg_or($spec, 'html_title', $project_tag . ' differential GWAS Manhattan Plot'),
        LOCAL_OUTPUT_PREFIX => $local_output_prefix,
        LOCAL_HTML_TITLE => cfg_or($spec, 'local_html_title', $project_tag . ' Local Top Hits Manhattan Plot'),
        LOCAL_TOP_HITS_CSV_BASENAME => $local_top_hits_csv_basename,
        OUTPUT_HTML_BASENAME => $output_html_basename,
        FOREST_OUTPUT_PREFIX => $forest_output_prefix,
        FOREST_HTML_TITLE => cfg_or($spec, 'forest_html_title', $project_tag . ' top-hit forest plots'),
        FOREST_TOP_HITS_CSV_BASENAME => $forest_top_hits_csv_basename,
        FOREST_OUTPUT_HTML_BASENAME => $forest_output_html_basename,
        FOREST_OUTPUT_MANIFEST_BASENAME => $forest_output_manifest_basename,
        PAIR_DEFS => $pair_info->{pair_defs},
        GROUP_TRACKS => $pair_info->{group_rep},
    };
}

sub suggest_genomewide_fig_height {
    my ($track_count) = @_;
    $track_count ||= 1;
    my $height = 260 + (110 * $track_count);
    $height = 420 if $height < 420;
    return int($height);
}

sub suggest_local_manhattan_fig_height {
    my ($track_count, $hits_per_fig) = @_;
    $track_count ||= 1;
    $hits_per_fig ||= 6;
    my $height = 260 + (115 * $track_count) + (55 * ($hits_per_fig > 6 ? 6 : $hits_per_fig));
    $height = 760 if $height < 760;
    return int($height);
}

sub build_alias_map {
    my ($pair_info) = @_;
    my %map = map { $_ . '_STD_P' => $_ . '_STD_DIFF_P' } @{ $pair_info->{prefix_order} };
    for my $group (sort keys %{ $pair_info->{group_rep} }) {
        my $label = safe_name($pair_info->{group_rep}{$group}{label});
        $map{$label . '_P'} = $pair_info->{group_rep}{$group}{pvar};
    }
    return \%map;
}

sub build_post_alias_map {
    my ($pair_info) = @_;
    my %map = map { $_ . '_STD_Z' => $_ . '_STD_DIFF_Z' } @{ $pair_info->{prefix_order} };
    for my $group (sort keys %{ $pair_info->{group_rep} }) {
        my $label = safe_name($pair_info->{group_rep}{$group}{label});
        $map{$label . '_Z'} = $pair_info->{group_rep}{$group}{zvar};
    }
    return \%map;
}

sub validate_spec {
    my ($spec) = @_;
    die "Spec root must be a JSON object\n" unless ref($spec) eq 'HASH';
    validate_alias_override_spec($spec->{raw_column_aliases}) if exists $spec->{raw_column_aliases};
    my $mode = cfg_or($spec, 'source_mode', 'raw_pgc_vcf_sumstats');
    die "Unsupported source_mode: $mode\n"
      unless $mode =~ /^(?:raw_pgc_vcf_sumstats|precomputed_diff|precomputed_diff_stdized|merged_gwas_table)$/;

    if ($mode eq 'raw_pgc_vcf_sumstats') {
        die "input_dir is required for raw_pgc_vcf_sumstats\n" unless length cfg_or($spec, 'input_dir', '');
        die "groups must be a non-empty array for raw_pgc_vcf_sumstats\n"
          unless ref($spec->{groups}) eq 'ARRAY' && @{ $spec->{groups} };
        my %group_tags;
        for my $g (@{ $spec->{groups} }) {
            die "Each group must be an object\n" unless ref($g) eq 'HASH';
            my $tag = $g->{tag} // '';
            die "Each group needs a tag\n" unless length $tag;
            die "Duplicate group tag: $tag\n" if $group_tags{$tag}++;
            die "Group $tag needs a non-empty files array\n"
              unless ref($g->{files}) eq 'ARRAY' && @{ $g->{files} };
        }
        for my $pair (@{ $spec->{pairs} || [] }) {
            die "Pair references unknown group1 $pair->{group1}\n" unless $group_tags{ $pair->{group1} // '' };
            die "Pair references unknown group2 $pair->{group2}\n" unless $group_tags{ $pair->{group2} // '' };
        }
    }

    if ($mode eq 'precomputed_diff') {
        die "input_diff is required for precomputed_diff\n" unless length cfg_or($spec, 'input_diff', '');
    }
    if ($mode eq 'precomputed_diff_stdized') {
        die "input_stdized is required for precomputed_diff_stdized\n" unless length cfg_or($spec, 'input_stdized', '');
    }
    if ($mode eq 'merged_gwas_table') {
        die "input_merged is required for merged_gwas_table\n" unless length cfg_or($spec, 'input_merged', '');
        die "merged_group_tracks must be a non-empty array for merged_gwas_table\n"
          unless ref($spec->{merged_group_tracks}) eq 'ARRAY' && @{ $spec->{merged_group_tracks} } >= 2;
    }
    die "pairs must be a non-empty array\n"
      unless ref($spec->{pairs}) eq 'ARRAY' && @{ $spec->{pairs} };
}

sub validate_alias_override_spec {
    my ($aliases) = @_;
    die "raw_column_aliases must be a JSON object\n" unless ref($aliases) eq 'HASH';
    for my $key (sort keys %{$aliases}) {
        my $val = $aliases->{$key};
        die "raw_column_aliases.$key must be a string or array of strings\n"
          if ref($val) && ref($val) ne 'ARRAY';
        my @vals = ref($val) eq 'ARRAY' ? @{$val} : ($val);
        die "raw_column_aliases.$key must not be empty\n" unless @vals;
        for my $alias (@vals) {
            die "raw_column_aliases.$key contains an empty alias\n"
              unless defined $alias && length $alias;
        }
    }
}

sub resolve_oda_helper_unix {
    my ($bin_dir) = @_;
    my @candidates = (
        File::Spec->catfile($bin_dir, 'DiffGWASDeps', 'run_sas_codes_or_script_in_ODA.pl'),
        File::Spec->catfile($bin_dir, 'run_sas_codes_or_script_in_ODA.pl'),
    );
    for my $candidate (@candidates) {
        next unless defined $candidate && -f $candidate;
        return normalize_unix_path($candidate);
    }
    return 'run_sas_codes_or_script_in_ODA.pl';
}

sub validate_generated_files {
    my ($generated, $pair_info) = @_;
    for my $path ($generated->{wide_output}, $generated->{wide_manifest}, $generated->{runner_config}, $generated->{preset_config}) {
        die "Expected generated file missing or empty: $path\n" unless defined $path && -s cygpath_to_win($path);
    }
    assert_wide_output_matches_pairs($generated->{wide_output}, $pair_info);
    assert_not_single_snp_manifest($generated->{wide_manifest});
    assert_wide_manifest_matches_pairs($generated->{wide_manifest}, $pair_info);
}

sub assert_not_single_snp_manifest {
    my ($manifest_path) = @_;
    #print STDERR "Manifest path is not defined\n" and;
    return unless defined $manifest_path && length $manifest_path;
    my $win = cygpath_to_win($manifest_path);
    #print STDERR "Manifest file does not exist or is empty: $win\n" and;
    return unless -e $win && -s $win;

    open my $fh, '<', $win or die "Cannot read manifest $manifest_path: $!\n";
    my %metric;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my ($k, $v) = split /\t/, $line, 2;
        next unless defined $k && defined $v;
        $metric{$k} = $v;
    }
    close $fh;

    if (exists $metric{target_snp} || exists $metric{rows_in_window}) {
        my $target = $metric{target_snp} // 'unknown_target';
        die "Wide subset manifest $manifest_path is a single-SNP/local-window artifact (target_snp=$target). Rebuild the full extract_wide_subset output before plotting or common-hit verification.\n";
    }
}

sub assert_wide_manifest_matches_pairs {
    my ($manifest_path, $pair_info) = @_;
    return unless defined $manifest_path && length $manifest_path;
    return unless ref($pair_info) eq 'HASH';
    my $win = cygpath_to_win($manifest_path);
    return unless -e $win && -s $win;

    open my $fh, '<', $win or die "Cannot read manifest $manifest_path: $!\n";
    my %metric;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my ($k, $v) = split /\t/, $line, 2;
        next unless defined $k && defined $v;
        $metric{$k} = $v;
    }
    close $fh;

    my @expected_prefixes = @{ $pair_info->{prefix_order} || [] };
    return unless @expected_prefixes;
    my %expected = map { $_ => 1 } @expected_prefixes;
    my %seen = map { $_ => 1 } grep { length } split /\s*,\s*/, ($metric{pair_prefixes} // '');

    my @missing_prefixes = grep { !$seen{$_} } @expected_prefixes;
    if (@missing_prefixes) {
        my $found = $metric{pair_prefixes} // '(none)';
        print STDERR "[warn] Wide subset manifest $manifest_path is stale or incomplete for the current spec. Missing pair prefix(es): "
          . join(', ', @missing_prefixes)
          . ". Found pair_prefixes=$found. The pipeline will trust the actual wide-output header for plotting.\n";
        return;
    }

    my %columns = map { $_ => 1 } grep { length } split /\s*,\s*/, ($metric{columns} // '');
    my @required_columns;
    for my $prefix (@expected_prefixes) {
        push @required_columns,
          "${prefix}_GROUP1_BETA",
          "${prefix}_GROUP2_BETA",
          "${prefix}_DIFF_BETA",
          "${prefix}_GROUP1_SE",
          "${prefix}_GROUP2_SE",
          "${prefix}_DIFF_SE",
          "${prefix}_GROUP1_P",
          "${prefix}_GROUP2_P",
          "${prefix}_DIFF_P",
          "${prefix}_STD_DIFF_Z",
          "${prefix}_STD_DIFF_P";
    }
    my @missing_columns = grep { !$columns{$_} } @required_columns;
    if (@missing_columns) {
        my $show = @missing_columns > 8 ? 8 : scalar @missing_columns;
        print STDERR "[warn] Wide subset manifest $manifest_path is missing required multi-pair columns for the current spec, for example: "
          . join(', ', @missing_columns[0 .. $show - 1])
          . ". The pipeline will trust the actual wide-output header for plotting.\n";
    }
}

sub assert_wide_output_matches_pairs {
    my ($wide_output, $pair_info) = @_;
    return unless defined $wide_output && length $wide_output;
    return unless ref($pair_info) eq 'HASH';

    my $win = cygpath_to_win($wide_output);
    die "Wide subset file missing or empty: $wide_output\n" unless -e $win && -s $win;

    my $wide_output_for_shell = $wide_output;
    if ($wide_output_for_shell !~ m{^(?:[A-Za-z]:)?[\\/]} && defined $workdir && length $workdir) {
        $wide_output_for_shell = normalize_unix_path(File::Spec->catfile($workdir, $wide_output_for_shell));
    }

    my @expected_prefixes = @{ $pair_info->{prefix_order} || [] };
    return unless @expected_prefixes;

    my @required_columns = qw(CHR BP A1 A2 SNP);
    for my $prefix (@expected_prefixes) {
        push @required_columns,
          "${prefix}_GROUP1_BETA",
          "${prefix}_GROUP2_BETA",
          "${prefix}_DIFF_BETA",
          "${prefix}_GROUP1_SE",
          "${prefix}_GROUP2_SE",
          "${prefix}_DIFF_SE",
          "${prefix}_GROUP1_P",
          "${prefix}_GROUP2_P",
          "${prefix}_DIFF_P",
          "${prefix}_STD_DIFF_Z",
          "${prefix}_STD_DIFF_P";
    }

    my $cmd = qq{"$bash_path" -lc 'gzip -cd "$wide_output_for_shell" | head -n 1'};
    my $header = qx{$cmd};
    die "Failed to read header from wide subset $wide_output\n" unless defined $header && length $header;
    chomp $header;
    $header =~ s/\r$//;
    my %columns = map { $_ => 1 } split /\t/, $header;
    my @missing = grep { !$columns{$_} } @required_columns;
    if (@missing) {
        my $show = @missing > 8 ? 8 : scalar @missing;
        die "Wide subset file $wide_output is missing required columns for the current spec, for example: "
          . join(', ', @missing[0 .. $show - 1])
          . ". Rerun extract_wide_subset with --force before plotting.\n";
    }
}

sub oda_remote_file_exists {
    my ($bash_path, $workdir, $session_id, $basename) = @_;
    return 0 unless defined $basename && length $basename;
    my $cmd = qq{"$bash_path" -lc 'cd "$workdir" && perl "$oda_helper_unix" --dir4listing "~" --output-prefix "check_remote_plot_data" 2>&1'};
    my $output = qx{$cmd};
    return ($output =~ /^\Q$basename\E$/m) ? 1 : 0;
}

sub write_json_if_defined {
    my ($path, $data) = @_;
    #print STDERR "Path is not defined\n" and ;
    return unless defined $path && defined $data;
    my $win = cygpath_to_win($path);
    ensure_parent_dir($win);
    open my $fh, '>', $win or die "Cannot write $win: $!\n";
    print {$fh} pretty_json($data), "\n";
    close $fh;
    #print STDERR "Wrote JSON to $path\n";
}

sub pretty_json {
    my ($data) = @_;
    my $json = JSON::PP->new->ascii->pretty->canonical->encode($data);
    return $json;
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

sub run_step {
    my (%args) = @_;
    my $name = $args{name};
    my $command = $args{command};
    my $outputs = $args{outputs} || [];
    my $force = $args{force} || 0;
    my $step_started = time();
    unless ($force) {
        my $all_exist = 1;
        for my $out (@{$outputs}) {
            my $win = cygpath_to_win($out);
            if (!defined $out || !length $out || !-s $win) {
                $all_exist = 0;
                last;
            }
        }
        if ($all_exist) {
            print "[skip] $name already has expected outputs\n";
            print "[done] $name finished in " . format_elapsed_seconds(time() - $step_started) . " (reused existing outputs)\n";
            return;
        }
    }
    print "[run] $name\n";
    my ($stderr_fh, $stderr_path) = tempfile('auto_prepare_step_stderr_XXXX', SUFFIX => '.log', TMPDIR => 1);
    close $stderr_fh;
    my $wrapped_command = qq{$command 2>"$stderr_path"};
    my $rc = system($wrapped_command);
    if ($rc != 0) {
        my $stderr_text = '';
        if (open my $err_fh, '<', $stderr_path) {
            local $/;
            $stderr_text = <$err_fh> // '';
            close $err_fh;
        }
        unlink $stderr_path if defined $stderr_path && -f $stderr_path;
        my $hint = latest_step_log_hint($name);
        my $msg = "Step $name failed: $command\n";
        if (defined $stderr_text && $stderr_text =~ /\S/) {
            $stderr_text =~ s/\s+\z//;
            $msg .= "Captured stderr:\n$stderr_text\n";
        }
        if (defined $hint && length $hint) {
            $msg .= "Latest step log hint: $hint\n";
        }
        die $msg;
    }
    unlink $stderr_path if defined $stderr_path && -f $stderr_path;
    print "[done] $name finished in " . format_elapsed_seconds(time() - $step_started) . "\n";
}

sub latest_step_log_hint {
    my ($step_name) = @_;
    my %patterns = (
        plot_local_manhattan => 'run_local_hits_manhattan_png_*',
        plot_local_gtf => 'run_local_hits_with_gtf_*',
        plot_manhattan => 'run_manhattan4diffgwas_png_*',
        plot_forest => 'run_top_hits_forest_plot_*',
    );
    my $pattern = $patterns{$step_name} || return '';
    my @dirs = grep { -d $_ } glob(File::Spec->catfile($Bin, $pattern));
    return '' unless @dirs;
    @dirs = sort { (stat($b))[9] <=> (stat($a))[9] } @dirs;
    my $log = File::Spec->catfile($dirs[0], 'output.html.info.txt');
    return -f $log ? $log : $dirs[0];
}

sub resolve_step_selection {
    my (%args) = @_;
    my @step_args = @{ $args{step_args} || [] };
    my %step_flag = %{ $args{step_flag} || {} };
    my $from_step = canonical_step_name($args{from_step} // '');
    my $to_step = canonical_step_name($args{to_step} // '');
    my @step_defs = @{ $args{step_defs} || [] };
    my @available = map { $_->{name} } @step_defs;
    my %available = map { $_ => 1 } @available;

    my @explicit_steps;
    for my $raw (@step_args) {
        push @explicit_steps, grep { length } map { canonical_step_name($_) } split /\s*,\s*/, ($raw // '');
    }
    for my $name (sort keys %step_flag) {
        push @explicit_steps, $name if $step_flag{$name};
    }
    my %explicit_requested = map { $_ => 1 } @explicit_steps;

    my %selected;
    my $selection_mode = 'full_pipeline';

    if (@explicit_steps) {
        $selection_mode = 'explicit_steps';
        for my $name (@explicit_steps) {
            die "Unknown step: $name\n" unless $available{$name};
            $selected{$name} = 1;
        }
    }
    elsif (length $from_step || length $to_step) {
        $selection_mode = 'step_range';
        die "Unknown from-step: $from_step\n" if length $from_step && !$available{$from_step};
        die "Unknown to-step: $to_step\n" if length $to_step && !$available{$to_step};
        my $seen_start = !length $from_step;
        for my $name (@available) {
            $seen_start = 1 if length $from_step && $name eq $from_step;
            next unless $seen_start;
            $selected{$name} = 1;
            last if length $to_step && $name eq $to_step;
        }
        if (length $from_step && !$selected{$from_step}) {
            die "from-step $from_step is not active for this source_mode/plot selection\n";
        }
        if (length $to_step && !$selected{$to_step}) {
            die "to-step $to_step is not active for this source_mode/plot selection\n";
        }
    }
    else {
        %selected = %available;
    }

    my @selected_order = grep { $selected{$_} } @available;
    return {
        selected       => \%selected,
        selected_order => \@selected_order,
        available      => \@available,
        selection_mode => $selection_mode,
        explicit_requested => \%explicit_requested,
    };
}

sub print_step_catalog {
    my ($step_defs, $step_selection) = @_;
    print "Available pipeline steps for the current spec:\n";
    for my $step (@{$step_defs}) {
        my $status = (exists $step->{enabled} && !$step->{enabled}) ? 'inactive' : 'active';
        print "  $step->{name} [$status]\n";
        print "    $step->{description}\n" if defined $step->{description} && length $step->{description};
    }
    print "\nCurrent selection mode: $step_selection->{selection_mode}\n";
    print "Selected steps in execution order:\n";
    print "  $_\n" for @{ $step_selection->{selected_order} };
}

sub canonical_step_name {
    my ($name) = @_;
    return '' unless defined $name;
    $name =~ s/^\s+|\s+$//g;
    $name =~ tr/A-Z/a-z/;
    $name =~ s/-/_/g;
    return $name;
}

sub print_summary {
    my ($summary) = @_;
    print "Automation complete.\n";
    for my $key (sort keys %{$summary}) {
        my $value = $summary->{$key};
        if (ref($value) eq 'ARRAY') {
            print "$key:\n";
            print "  $_\n" for @{$value};
        }
        else {
            print "$key: $value\n";
        }
    }
}

sub ensure_parent_dir {
    my ($path) = @_;
    my ($vol, $dir) = File::Spec->splitpath($path);
    my $parent = File::Spec->catpath($vol, $dir, '');
    return unless length $parent;
    mkdir $parent unless -d $parent;
}

sub load_json {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    local $/;
    my $json = <$fh>;
    close $fh;
    my $cfg = decode_json($json);
    die "JSON root must be an object: $path\n" unless ref($cfg) eq 'HASH';
    return $cfg;
}

sub load_alias_override_file {
    my ($path) = @_;
    return {} unless defined $path && length $path;
    my $cfg = load_json($path);
    validate_alias_override_spec($cfg);
    return $cfg;
}

sub merge_alias_override_specs {
    my ($base, $extra) = @_;
    $base = {} unless ref($base) eq 'HASH';
    $extra = {} unless ref($extra) eq 'HASH';
    my %merged = %{$base};
    for my $key (keys %{$extra}) {
        my @vals = ref($extra->{$key}) eq 'ARRAY' ? @{ $extra->{$key} } : ($extra->{$key});
        my @base_vals = exists $merged{$key}
          ? (ref($merged{$key}) eq 'ARRAY' ? @{ $merged{$key} } : ($merged{$key}))
          : ();
        my %seen;
        $merged{$key} = [
            grep { defined $_ && length $_ && !$seen{ normalize_header_name($_) }++ }
            (@vals, @base_vals)
        ];
    }
    return \%merged;
}

sub cfg_or {
    my ($cfg, $key, $fallback) = @_;
    return $fallback unless exists $cfg->{$key};
    return $cfg->{$key};
}

sub script_root_dir {
    my $dir = eval { abs_path($Bin) } || File::Spec->rel2abs($Bin) || $Bin;
    return normalize_unix_path($dir);
}

sub resolve_pipeline_workdir {
    my (%args) = @_;
    my $spec = $args{spec} || {};
    my $fallback = normalize_unix_path($args{fallback} // $Bin);

    if (defined $ENV{PIPELINE_WORKDIR} && length $ENV{PIPELINE_WORKDIR}) {
        my $env_workdir = normalize_unix_path($ENV{PIPELINE_WORKDIR});
        $spec->{workdir} = $env_workdir if ref($spec) eq 'HASH';
        return $env_workdir;
    }

    my $configured = normalize_unix_path(cfg_or($spec, 'workdir', $fallback));
    return $configured if $ENV{PIPELINE_RESPECT_SPEC_WORKDIR};

    if (
        length($configured)
        && length($fallback)
        && $configured ne $fallback
      )
    {
        warn "[info] Overriding spec workdir '$configured' with current pipeline directory '$fallback' for portable execution.\n";
        $spec->{workdir} = $fallback if ref($spec) eq 'HASH';
        return $fallback;
    }

    $spec->{workdir} = $configured if ref($spec) eq 'HASH';
    return $configured;
}

sub normalize_target_snp_gene_overrides {
    my ($raw) = @_;
    return '' unless defined $raw;
    if (ref($raw) eq 'HASH') {
        return join(
            ',',
            map {
                my $gene = $raw->{$_};
                defined($gene) && length($gene)
                  ? ($_.':' . $gene)
                  : ()
            } sort keys %{$raw}
        );
    }
    return $raw if !ref($raw);
    return '';
}

sub population_map_to_string {
    my ($raw) = @_;
    return '' unless defined $raw;
    if (ref($raw) eq 'HASH') {
        return join(
            ',',
            map {
                my $value = $raw->{$_};
                defined($value) && length($value) ? ($_. '=' . $value) : ()
            } sort keys %{$raw}
        );
    }
    return $raw if !ref($raw);
    return '';
}

sub safe_name {
    my ($text) = @_;
    $text = 'diff_gwas' unless defined $text && length $text;
    $text =~ s/[^\w.+-]+/_/g;
    $text =~ s/_+/_/g;
    $text =~ s/^_+|_+$//g;
    return $text;
}

sub trim {
    my ($value) = @_;
    return '' unless defined $value;
    $value =~ s/^\s+//;
    $value =~ s/\s+$//;
    return $value;
}

sub differential_threshold_fallback_value {
    my ($primary, $fallback) = @_;
    $primary = trim($primary);
    $fallback = trim($fallback);
    return $fallback if length $fallback;
    return '1e-5' if normalized_numeric_text($primary) eq normalized_numeric_text('1e-6');
    return '';
}

sub build_differential_threshold_ladder {
    my ($primary, $fallback) = @_;
    $primary = trim($primary);
    my $resolved_fallback = differential_threshold_fallback_value($primary, $fallback);
    my @vals;
    for my $value ($primary, $resolved_fallback) {
        next unless defined $value && length $value;
        push @vals, $value unless grep { $_ eq $value } @vals;
    }
    return join(' ', @vals);
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

sub infer_group_label {
    my ($tag) = @_;
    return '' unless defined $tag;
    my $label = $tag;
    $label =~ s/^.*?_W\d+_//;
    $label =~ s/^SCZ_//;
    $label =~ s/_/ /g;
    return $label;
}

sub normalize_unix_path {
    my ($path) = @_;
    return '' unless defined $path;
    my $is_unc = ($path =~ m{^(?:\\\\|//)[^\\/]+[\\/]} ) ? 1 : 0;
    $path =~ s{\\}{/}g;
    $path =~ s{//+}{/}g;
    $path = '/' . $path if $is_unc && $path !~ m{^//};
    $path = '/' . $path if $is_unc && $path !~ m{^//};
    return $path;
}

sub cygpath_to_win {
    my ($path) = @_;
    #print STDERR "Path is not defined\n" and ;
    return '' unless defined $path;
    #print STDERR "Path is defined as $path\n" and ;
    return $path  if $path =~ /^[A-Za-z]:[\\\/]/;
    return $path if $^O !~ /^(?:cygwin|MSWin32)$/i;
    if ($path =~ m{^/mnt/([A-Za-z])/(.*)$}) {
        my ($drive, $rest) = ($1, $2);
        $rest =~ s{/}{\\}g;
        #print STDERR "Path is defined as $path and converted to $drive:\\$rest\n" and; 
        return uc($drive) . ":\\" . $rest;
    }
    if ($path =~ m{^/cygdrive/([A-Za-z])/(.*)$}) {
        my ($drive, $rest) = ($1, $2);
        $rest =~ s{/}{\\}g;
        return uc($drive) . ":\\" . $rest;
    }
    my $win = $path;
    $win =~ s{/}{\\}g;
    #print STDERR "Path is defined as $path and converted to $win\n" and ;
    return $win;
}

sub usage {
    return <<"USAGE";
Usage:
  perl auto_prepare_and_run_diff_gwas.pl --spec comparison_spec.json [options]
  perl auto_prepare_and_run_diff_gwas.pl --gwas-dir /mnt/e/path/to/gwas_dir --generate-spec-only [options]

Options:
  --spec FILE.json      Required comparison spec
  --gwas-dir DIR        Scan a directory of potential raw GWAS summary-statistics files and draft a spec JSON
  --spec-out FILE.json  Where to write the auto-generated spec JSON
  --raw-column-alias-config FILE.json
                       Optional JSON mapping canonical raw columns like ID/POS/PVAL
                       to extra header aliases for new GWAS formats
  --preview-spec        Print the inferred spec JSON to STDOUT before writing it
  --mode MODE           full|configs . Default: full
  --skip-plots          Generate/prepare data but do not call SAS ODA runners
  --plots LIST          Comma-separated plot set. Default: manhattan,local_manhattan,local_gtf
  --force               Rerun steps even if expected outputs already exist
  --list-steps          Print the available step names for the current spec and exit
  --step NAME           Run only one named step; repeat or comma-separate to run several
  --from-step NAME      Run from a named step onward
  --to-step NAME        Stop after a named step
  --merge-raw           Convenience alias for --step merge_raw
  --sort-long           Convenience alias for --step sort_long
  --diff-pairs          Convenience alias for --step diff_pairs
  --standardize-diff    Convenience alias for --step standardize_diff
  --extract-wide-subset Convenience alias for --step extract_wide_subset
  --plot-manhattan      Convenience alias for --step plot_manhattan
  --plot-local-manhattan Convenience alias for --step plot_local_manhattan
  --plot-local-gtf      Convenience alias for --step plot_local_gtf
  --plot-forest         Convenience alias for --step plot_forest
  --emit-local-sas-scripts
                      Also emit local desktop-SAS runnable plot scripts for
                      plot_manhattan, plot_local_manhattan, plot_local_gtf,
                      and plot_forest.
  --local-sas-only     Emit the local desktop-SAS plot scripts and stop before
                      any SAS ODA submit/upload work for plot steps.
  --local-max-hits-per-fig N
                       Requested upper bound for local top-hit columns per panel.
                       Current pipeline maximum is 15.
  --local-manhattan-angle4xaxis-label N
                       Override the SNP/gene label rotation angle in local
                       top-hit Manhattan panels. Default is macro-driven.
  --local-manhattan-xgrp-y-pos N
                       Override the vertical position of SNP/gene labels in
                       local top-hit Manhattan panels.
  --local-manhattan-yoffset-top N
                       Override the top y-axis offset used to make room for
                       local Manhattan SNP/gene labels.
  --local-manhattan-yoffset-bottom N
                       Override the bottom y-axis offset paired with the top
                       offset for local Manhattan label layout.
  --local-manhattan-fontsize N
                       Override the base SAS font size used by the local
                       Manhattan macro, including SNP/gene label sizing.
  --local-manhattan-y-axis-label-size N
                       Override the y-axis title font size for local Manhattan.
  --local-manhattan-y-axis-value-size N
                       Override the y-axis tick font size for local Manhattan.
  --target-snps rs1,rs2,rs3
  --target-snp-genes rs1:GENE1,rs2:GENE2
                       Use an explicit SNP list, in the provided order, for the
                       local Manhattan and local GTF plot stages instead of
                       automatic top-hit picking.
  --display-gwas LIST
                       Optional comma-separated GWAS track selection shared by
                       both SAS ODA and gunplot pipelines. Use pair prefixes
                       such as ALL,EUR,ASN for differential tracks and GWAS
                       labels such as ALL_FEMALE or EUR_MALE for single-GWAS
                       tracks. When only one GWAS track is selected, the
                       pipeline renders single-GWAS genomewide/local Manhattan
                       and local GTF plots and uses that selected GWAS for
                       top-hit selection unless --target-snps is provided.
  --reference-build BUILD
                       Optional genome-build override for local gene-track
                       annotations. Accepted values: hg19, hg38, t2t. When
                       omitted, the pipeline first looks for build tokens in
                       input filenames and headers, then falls back to hg38.
  --sas-oda-account EMAIL
  --sas-oda-password PASS
                       Optional noninteractive SAS ODA credential bootstrap for
                       the first SAS-backed run. If omitted, the vendored SAS
                       helper now prompts interactively on first use, validates
                       the login with PROC SETINIT, and saves the working
                       account/password to the SASPy authinfo file.
  --prompt-sas-oda-auth
                       Force a SAS ODA credential refresh on the next helper
                       connection, even if an authinfo entry already exists.
  --local-gtf-label-snps rs1,rs2,rs3
                       Label these SNP names on top of the local GTF figure.
                       Default behavior labels the target/top-hit SNPs.
  --local-gtf-label-layout MODE
                       local GTF top-label layout: auto|vertical|horizontal.
  --local-gtf-window-bp BP
                       Override the genomic half-window used only for the local
                       GTF plot stage. This does not change local Manhattan.
  --local-gtf-yaxis-offset4max N
                       Set the starting local-GTF top headroom fraction passed
                       to the lattice macro as yaxis_offset4max. This matters
                       most for a single top SNP label and may still be updated
                       internally by the pipeline's SAS auto-tuning logic.
  --local-gtf-yoffset4textlabels N
                       Set the starting local-GTF top SNP label shift passed to
                       the lattice macro as Yoffset4textlabels. This moves the
                       SNP label up/down in y-axis-value units and may still be
                       adjusted internally by the pipeline.
  --exclude-non-protein-coding-genes-in-local-gtf
                       Legacy convenience flag. Local GTF bottom gene tracks are
                       protein-coding-only by default. To include non-coding
                       genes, set include_non_protein_coding_genes_in_local_gtf
                       to 1 in the spec JSON.
  --get-common-associations[=THR]
                       Select local top hits from replicable single-GWAS
                       associations that also show nominal association in
                       another GWAS with the same effect direction.
                       Optional THR sets the starting common-association
                       top-hit threshold for any single-GWAS association P;
                       default ladder starts at 5e-8.
  --common-association-top-hit-threshold THR
                       Explicit alias for the starting common-association
                       top-hit threshold for any single-GWAS association P.
  --cleanup-shared-plot-data Convenience alias for --step cleanup_shared_plot_data
  --generate-spec-only  When used with --gwas-dir, write or preview the draft spec and exit
  --project-tag TEXT    Optional override for inferred project_tag during auto-spec generation
  --artifact-stem TEXT  Optional override for inferred artifact_stem during auto-spec generation
  --print-spec-example  Print a full JSON spec example and exit
  --print-columns-help  Print required-column help and exit

The spec can describe:
  - raw_pgc_vcf_sumstats    : merge raw group files, sort, diff, standardize, extract, plot
  - precomputed_diff        : standardize, extract, plot
  - precomputed_diff_stdized: extract, plot
  - merged_gwas_table       : normalize one merged-wide GWAS table, then plot

Recommended interactive pattern:
  - open this repository in VS Code
  - on Windows, use the portable Cygwin integrated terminal
  - on macOS or Ubuntu Linux, use the native integrated shell
  - keep Codex and the local shell in the same workspace when driving the
    pipeline through the MCP server
USAGE
}

sub full_help {
    my (%args) = @_;
    my $show_example = exists $args{show_example} ? $args{show_example} : 1;
    my $show_columns = exists $args{show_columns} ? $args{show_columns} : 1;

    my $text = usage();
    $text .= "\nQuick start:\n";
    $text .= "  1. Prepare a JSON spec file describing your comparison.\n";
    $text .= "  2. Run:\n";
    $text .= "       perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json\n";
    $text .= "  3. Or only generate the auto-configs first:\n";
    $text .= "       perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --mode configs\n";
    $text .= "  4. Or let the script draft a spec from a GWAS directory:\n";
    $text .= "       perl auto_prepare_and_run_diff_gwas.pl --gwas-dir /mnt/e/path/to/gwas_dir --generate-spec-only\n";
    $text .= "  5. Or preview the inferred spec JSON without writing a file:\n";
    $text .= "       perl auto_prepare_and_run_diff_gwas.pl --gwas-dir /mnt/e/path/to/gwas_dir --preview-spec --generate-spec-only\n";
    $text .= "  6. If header auto-detection misses a new GWAS format, provide alias overrides:\n";
    $text .= "       perl auto_prepare_and_run_diff_gwas.pl --gwas-dir /mnt/e/path/to/gwas_dir --raw-column-alias-config ./raw_aliases.json --generate-spec-only\n";
    $text .= "\nEditor-centered usage:\n";
    $text .= "  - the same command works well from a VS Code workspace where Codex is active\n";
    $text .= "  - on Windows, prefer the portable Cygwin terminal opened from the repository\n";
    $text .= "  - on macOS and Ubuntu Linux, use the native integrated shell in the same workspace\n";
    $text .= "  - when using the MCP route, keep server.pl running in one integrated terminal and\n";
    $text .= "    let Codex orchestrate the same local scripts from that workspace\n";
    $text .= "\nWhat the script does in full mode:\n";
    $text .= "  - validates the comparison spec\n";
    $text .= "  - writes auto-generated merge/diff/preset/runner JSON configs\n";
    $text .= "  - reuses existing outputs unless --force is supplied\n";
    $text .= "  - prepares differential GWAS artifacts\n";
    $text .= "  - runs the requested SAS ODA plotting wrappers\n";
    $text .= "  - on the first SAS ODA run, the helper now prompts for the SAS ODA\n";
    $text .= "    account/password if they are not already saved, validates them with\n";
    $text .= "    proc setinit;run;, and then reuses the saved authinfo entry later\n";
    $text .= "\nStep-aware reruns:\n";
    $text .= "  - use --list-steps to see the exact pipeline step names for the current spec\n";
    $text .= "  - use --step plot_local_gtf when you only want one stage rerun\n";
    $text .= "  - use --step plot_forest when you want OR forest plots for the current\n";
    $text .= "    top-hit list across the selected single-GWAS association tracks\n";
    $text .= "  - use --from-step extract_wide_subset to rerun from a middle stage onward\n";
    $text .= "  - the convenience flags like --plot-manhattan and --extract-wide-subset map to the same step names\n";
    $text .= "  - use --local-max-hits-per-fig N when you want narrower local-Manhattan panels;\n";
    $text .= "    the CLI/spec layer accepts up to 10 columns per panel before dense-track auto-tightening\n";
    $text .= "  - use --local-gtf-window-bp BP when you want the local GTF plot to use a\n";
    $text .= "    custom genomic half-window without changing the local Manhattan stage\n";
        $text .= "  - use --target-snps rs1,rs2,... when you want the local Manhattan and\n";
        $text .= "  - use --target-snp-genes rs1:GENE1,rs2:GENE2 to override adjacent-gene labels\n";
    $text .= "    local GTF plots to render an explicit SNP list in your chosen order\n";
    $text .= "  - use --local-gtf-label-snps rs1,rs2,... when you want to label extra or\n";
    $text .= "    specific SNP names on top of the local GTF plot window\n";
    $text .= "  - use --local-gtf-label-layout vertical or --local-gtf-label-layout horizontal\n";
    $text .= "    to choose how those top SNP labels are drawn; auto keeps the macro-driven default\n";
    $text .= "  - use --local-gtf-yaxis-offset4max N when you want to set the starting\n";
    $text .= "    top headroom fraction for local GTF plots; this is the parameter that\n";
    $text .= "    controls single-SNP top headroom after the macro overrides\n";
    $text .= "    yoffset4max_drawmarkersontop, although the pipeline may still auto-tune it\n";
    $text .= "  - use --local-gtf-yoffset4textlabels N when you want to set the starting\n";
    $text .= "    top SNP label shift in y-axis units for local GTF plots; the pipeline may\n";
    $text .= "    still adjust it internally for readability\n";
    $text .= "  - local GTF bottom tracks are protein-coding-only by default; set\n";
    $text .= "    include_non_protein_coding_genes_in_local_gtf=1 in the spec when you\n";
    $text .= "    explicitly want non-coding genes included\n";
    $text .= "  - use --get-common-associations when you want local top hits driven by\n";
    $text .= "    strong single-GWAS association first, then retained only when another GWAS\n";
    $text .= "    shows nominal association with the same effect direction\n";
    $text .= "  - you can optionally pass a starting threshold, for example\n";
    $text .= "    --get-common-associations=5e-8; this applies to single-GWAS association P,\n";
    $text .= "    and the ladder still relaxes automatically if needed\n";
    $text .= "\nPerformance note:\n";
    $text .= "  - set keep_remote_plot_data=1 in the spec when you expect to rerun plots\n";
    $text .= "    for the same comparison. That keeps the wide gz subset in SAS ODA and lets\n";
    $text .= "    later runs skip both re-upload and delete of the large plot input file.\n";

    if ($show_columns) {
        $text .= "\nRequired columns by source_mode:\n";
        $text .= "  raw_pgc_vcf_sumstats\n";
        $text .= "    Expected raw GWAS columns in each input file:\n";
        $text .= "      CHROM  chromosome\n";
        $text .= "      ID     SNP identifier, usually rsID\n";
        $text .= "      POS    base-pair position\n";
        $text .= "      A1     effect allele\n";
        $text .= "      A2     non-effect/reference allele as supplied by the source\n";
        $text .= "      FCAS   case/effect-allele frequency\n";
        $text .= "      FCON   control/effect-allele frequency\n";
        $text .= "      IMPINFO imputation/info score\n";
        $text .= "      BETA   GWAS effect estimate\n";
        $text .= "      SE     standard error of BETA\n";
        $text .= "      PVAL   association P-value\n";
        $text .= "    Optional raw GWAS columns currently tolerated:\n";
        $text .= "      NCAS NCON NEFF\n";
        $text .= "    Notes:\n";
        $text .= "      - this mode is designed for the PGC VCF-like summary-statistics tables\n";
        $text .= "      - all input GWAS files in one comparison should use the same coordinate and allele convention\n";
        $text .= "      - groups are defined in the spec under groups[].files\n";
        $text .= "      - when using --gwas-dir, the script scans candidate files, checks the header row,\n";
        $text .= "        and only includes files where all required fields can be resolved\n";
        $text .= "      - if there are not enough raw GWAS groups, the script will also try to detect\n";
        $text .= "        precomputed differential tables and draft a precomputed_diff or\n";
        $text .= "        precomputed_diff_stdized spec instead\n";
        $text .= "      - the raw merge helper will try to auto-detect common header synonyms\n";
        $text .= "        for required fields, for example:\n";
        $text .= "          chromosome: CHROM, CHR, CHROMOSOME\n";
        $text .= "          SNP id    : ID, SNP, RSID, RS, MARKERNAME\n";
        $text .= "          position  : POS, BP, POSITION, PS\n";
        $text .= "          P value   : PVAL, P, PVALUE, P_VALUE, p-value, P_LRT, P_WALD\n";
        $text .= "          A1        : A1, EA, EFFECT_ALLELE, ALLELE1\n";
        $text .= "          A2        : A2, NEA, OTHER_ALLELE, ALLELE2, REF, ALLELE0\n";
        $text .= "      - GEMMA association outputs with headers like chr/rs/ps/allele1/allele0/beta/se/p_lrt are accepted\n";
        $text .= "      - when GEMMA provides p_wald, p_lrt, and p_score together, p_lrt is used by default\n";
        $text .= "      - case/control frequencies, info, and sample-size columns are optional\n";
        $text .= "      - if a file uses unfamiliar header names, you can provide raw_column_aliases in the spec\n";
        $text .= "        or pass --raw-column-alias-config aliases.json when generating or rerunning the spec\n";
        $text .= "      - alias JSON format example:\n";
        $text .= "          {\"ID\":[\"MARKER_ID\"],\"POS\":[\"GENOMIC_POS\"],\"PVAL\":[\"PVALUE_LRT\"]}\n";
        $text .= "\n";
        $text .= "  precomputed_diff\n";
        $text .= "    Expected columns in input_diff:\n";
        $text .= "      CHR BP A1 A2 SNP PAIR_TAG\n";
        $text .= "      GROUP1_BETA GROUP2_BETA DIFF_BETA\n";
        $text .= "      GROUP1_SE   GROUP2_SE   DIFF_SE\n";
        $text .= "      GROUP1_P    GROUP2_P    DIFF_P\n";
        $text .= "      DIFF_Z\n";
        $text .= "    Optional metadata columns are fine and will be preserved downstream only when selected by configs.\n";
        $text .= "\n";
        $text .= "  precomputed_diff_stdized\n";
        $text .= "    Expected columns in input_stdized:\n";
        $text .= "      all columns from precomputed_diff plus:\n";
        $text .= "      STD_DIFF_Z STD_DIFF_P\n";
        $text .= "    This is the easiest mode when you already have a finalized standardized long differential table.\n";
        $text .= "\n";
        $text .= "  merged_gwas_table\n";
        $text .= "    Expected columns in input_merged:\n";
        $text .= "      shared locus columns such as CHR BP SNP A1 A2\n";
        $text .= "      at least two cohort-level BETA/SE/P blocks, for example:\n";
        $text .= "        BETA_DS_ALL SE_DS_ALL P_DS_ALL\n";
        $text .= "        BETA_MP2PRT SE_MP2PRT P_MP2PRT\n";
        $text .= "      optional extra association tracks with P/Z columns, for example:\n";
        $text .= "        PR_meta WEIGHTED_Z_meta\n";
        $text .= "    This mode derives pairwise differential columns automatically and can\n";
        $text .= "    also keep extra association tracks like meta-analysis P/Z columns for plotting.\n";
        $text .= "\n";
        $text .= "Meaning of the key differential columns:\n";
        $text .= "  GROUP1_BETA / GROUP2_BETA : effect sizes from the two GWASs being compared\n";
        $text .= "  DIFF_BETA                 : GROUP1_BETA - GROUP2_BETA\n";
        $text .= "  DIFF_SE                   : SE of the differential effect size\n";
        $text .= "  DIFF_Z                    : raw differential Z-score\n";
        $text .= "  DIFF_P                    : raw differential P-value\n";
        $text .= "  STD_DIFF_Z                : standardized differential Z-score used for cross-comparison plotting\n";
        $text .= "  STD_DIFF_P                : standardized differential P-value used for default plotting focus\n";
    }

    if ($show_example) {
        $text .= "\nExample comparison spec for ancestry-style pairwise GWAS comparisons:\n\n";
        $text .= sample_spec_json();
        $text .= "\n";
        $text .= "Field notes for the spec:\n";
        $text .= "  project_tag\n";
        $text .= "    Short tag used for plot output names.\n";
        $text .= "  artifact_stem\n";
        $text .= "    Base name used for intermediate/output tables on disk.\n";
        $text .= "  groups[].tag\n";
        $text .= "    Internal GWAS name used in pair definitions, for example SCZ_W3_ASN.\n";
        $text .= "  groups[].files\n";
        $text .= "    One or more raw files that belong to that GWAS group.\n";
        $text .= "  pairs[].pair_tag\n";
        $text .= "    Name written into the long differential table, for example SCZ_W3_ASN_vs_EUR.\n";
        $text .= "  pairs[].group1 / pairs[].group2\n";
        $text .= "    The two GWAS group tags being contrasted.\n";
        $text .= "  pairs[].prefix\n";
        $text .= "    Short prefix used in the wide table and SAS plotting variables, for example ASN_EUR.\n";
        $text .= "  pairs[].label\n";
        $text .= "    Human-readable label used in plot track names.\n";
        $text .= "  top_hit_focus_prefix\n";
        $text .= "    Which pairwise comparison should drive the local top-hit selection by default.\n";
        $text .= "  reference_build\n";
        $text .= "    Optional genome-build hint for local gene-track annotations. Accepted values are hg19,\n";
        $text .= "    hg38, and t2t. If omitted, the pipeline heuristically checks input filenames and headers,\n";
        $text .= "    then falls back to hg38. Explicitly set this field when your files do not carry a clear build token.\n";
        $text .= "  keep_remote_plot_data\n";
        $text .= "    When 1, keep the uploaded wide gz subset in SAS ODA for reuse across later plot reruns.\n";
        $text .= "  include_non_protein_coding_genes_in_local_gtf\n";
        $text .= "    Default is 0, so local GTF bottom tracks stay protein-coding-only. Set to 1 only when you explicitly want non-coding genes included.\n";
        $text .= "  local_manhattan_angle4xaxis_label / local_manhattan_xgrp_y_pos\n";
        $text .= "    Optional manual overrides for local Manhattan SNP/gene label rotation and vertical position.\n";
        $text .= "  local_manhattan_yoffset_top / local_manhattan_yoffset_bottom\n";
        $text .= "    Optional manual y-axis offsets used to make room for local Manhattan SNP/gene labels.\n";
        $text .= "  local_manhattan_fontsize / local_manhattan_y_axis_label_size / local_manhattan_y_axis_value_size\n";
        $text .= "    Optional manual text-size overrides for difficult local Manhattan label layouts.\n";
        $text .= "  gtf_yaxis_offset4max\n";
        $text .= "    Optional starting top headroom fraction for local GTF plots. In the patched\n";
        $text .= "    lattice macro this is the key single-SNP headroom control and can override\n";
        $text .= "    gtf_yoffset4max_drawmarkersontop after internal SAS auto-tuning.\n";
        $text .= "  gtf_yoffset4textlabels\n";
        $text .= "    Optional starting local-GTF SNP-label shift in y-axis units. This moves the\n";
        $text .= "    label itself up/down and can still be adjusted internally by the pipeline.\n";
        $text .= "\n";
        $text .= "Common demo commands:\n";
        $text .= "  Draft a spec from a directory and stop:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --gwas-dir /mnt/e/path/to/gwas_dir --generate-spec-only\n";
        $text .= "  Preview the inferred spec JSON without writing a file:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --gwas-dir /mnt/e/path/to/gwas_dir --preview-spec --generate-spec-only\n";
        $text .= "  Draft a spec to a chosen path:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --gwas-dir /mnt/e/path/to/gwas_dir --spec-out ./configs/my_auto_spec.json --generate-spec-only\n";
        $text .= "  Generate configs only:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --mode configs\n";
        $text .= "  Full pipeline:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json\n";
        $text .= "  Only rerun selected plots:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --plots manhattan,local_gtf\n";
        $text .= "  List the exact step names available for one spec:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --list-steps\n";
        $text .= "  Rerun only one stage:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_local_gtf\n";
        $text .= "  Note: local GTF now defaults to one locus per figure or batch unless you explicitly raise it with local_gtf_max_hits_per_fig in the spec or --local-max-hits-per-fig on the CLI.\n";
        $text .= "  Rerun one local GTF plot with manual top-label starting values:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_local_gtf --target-snps rs123 --local-gtf-yaxis-offset4max 0.02 --local-gtf-yoffset4textlabels 2\n";
        $text .= "  Rerun from a middle stage onward:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --from-step extract_wide_subset\n";
        $text .= "  Request narrower local top-hit panels:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_local_manhattan --local-max-hits-per-fig 4\n";
        $text .= "  Allow up to 15 local Manhattan columns in one panel:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_local_manhattan --local-max-hits-per-fig 15\n";
        $text .= "  Manually shift local Manhattan SNP/gene labels upward:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_local_manhattan --local-manhattan-xgrp-y-pos -3.4 --local-manhattan-yoffset-top 14 --local-manhattan-yoffset-bottom 0.5\n";
        $text .= "  Rotate local Manhattan SNP/gene labels and enlarge their base font:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_local_manhattan --local-manhattan-angle4xaxis-label 60 --local-manhattan-fontsize 3.0\n";
        $text .= "  Use a custom genomic window only for the local GTF plot:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_local_gtf --local-gtf-window-bp 2e7\n";
        $text .= "  Plot an explicit SNP list instead of auto-picked top hits:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_local_gtf --target-snps rs123,rs456,rs789\n";
        $text .= "  Label multiple SNPs on top of one local GTF window with vertical text:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_local_gtf --target-snps rs123 --local-gtf-label-snps rs123,rs456,rs789 --local-gtf-label-layout vertical\n";
        $text .= "  Label multiple SNPs on top of one local GTF window with horizontal text:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_local_gtf --target-snps rs123 --local-gtf-label-snps rs123,rs456,rs789 --local-gtf-label-layout horizontal\n";
        $text .= "  Explicitly include non-protein-coding genes in the local GTF bottom track by editing the spec:\n";
        $text .= "    \"include_non_protein_coding_genes_in_local_gtf\": 1\n";
        $text .= "  Use common/replicable association mode for local hits:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_local_manhattan --get-common-associations\n";
        $text .= "  Use common/replicable association mode with a custom starting threshold:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_local_manhattan --get-common-associations=5e-8\n";
        $text .= "  Render a focused forest plot for the top differential locus:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_forest --top-hit-max-loci 1\n";
        $text .= "  Render a forest plot for explicit targets in the pooled female/male GWASs only:\n";
        $text .= "    perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_forest --target-snps rs185665940 --display-gwas ALL_FEMALE,ALL_MALE\n";
    }

    return $text;
}

sub sample_spec_json {
    return <<'JSON';
{
  "source_mode": "raw_pgc_vcf_sumstats",
  "project_tag": "PGC_SCZ_ancestry_diff",
  "artifact_stem": "PGC_SCZ_ancestry_diff_effects",
  "input_dir": "/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Ancestry_Stratified_GWASs",
  "output_dir": "/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Ancestry_Stratified_GWASs",
  "workdir": "/mnt/g/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper",
  "reference_build": "hg19",
  "threshold": 0.05,
  "keep_remote_plot_data": 1,
  "top_hit_focus_prefix": "ASN_EUR",
  "top_hit_signal_thrshd": "1e-6",
  "local_window_bp": "1e7",
  "local_gtf_window_bp": "1e7",
  "include_non_protein_coding_genes_in_local_gtf": 0,
  "groups": [
    {
      "tag": "SCZ_W3_AFR",
      "files": ["PGC3_SCZ_wave3.afram.autosome.public.v3.vcf.tsv.gz"]
    },
    {
      "tag": "SCZ_W3_ASN",
      "files": ["PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.gz"]
    },
    {
      "tag": "SCZ_W3_EUR",
      "files": ["PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"]
    }
  ],
  "pairs": [
    {
      "pair_tag": "SCZ_W3_ASN_vs_EUR",
      "group1": "SCZ_W3_ASN",
      "group2": "SCZ_W3_EUR",
      "prefix": "ASN_EUR",
      "label": "ASN vs EUR"
    },
    {
      "pair_tag": "SCZ_W3_AFR_vs_EUR",
      "group1": "SCZ_W3_AFR",
      "group2": "SCZ_W3_EUR",
      "prefix": "AFR_EUR",
      "label": "AFR vs EUR"
    }
  ]
}
JSON
}
