#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use JSON::PP qw(decode_json encode_json);
use Test::More;
use lib "$Bin/../DiffGWASDeps";
use StudyDesignFilters qw(cross_ancestry_filter_advisory);

my $dir = tempdir('multigwas anemia XXXXX', TMPDIR => 1, CLEANUP => 1);
for my $name ('SAA_EAS', 'SAA_EUR_plus_EAS') {
    open my $out, '>', "$dir/$name.assoc.txt" or die $!;
    print {$out} "SNP\tchr\tpos\talt\tref\tbeta\tse\tp\n";
    print {$out} "rsPalindromic\t3\t108377251\tC\tG\t0.3\t0.2\t0.15\n";
    print {$out} "rsOrdinary\t3\t108377350\tA\tC\t0.4\t0.1\t0.0001\n";
    close $out or die $!;
}

my $driver = "$Bin/../auto_prepare_and_run_diff_gwas.pl";
my $default_spec_path = "$dir/anemia.default.spec.json";
is(system($^X, $driver, '--gwas-dir', $dir, '--spec-out', $default_spec_path,
    '--generate-spec-only'), 0, 'anemia inputs generate a spec without overrides');
my $default_spec = read_json($default_spec_path);
is($default_spec->{compute_meta}, 0, 'computed meta-analysis is off by default');
is($default_spec->{max_eaf_abs_diff}, '', 'frequency-discordance filtering is off by default');
like(join(' ', @{ $default_spec->{auto_detect_summary}{warnings} || [] }),
    qr/Cross-ancestry comparison.*max_eaf_abs_diff=none.*top_hit_maf_threshold=0.*population-specific/i,
    'inferred cross-ancestry spec warns about both default filters and real population-specific signals');
$default_spec->{configs_dir} = "$dir/default_configs";
open my $default_out, '>', $default_spec_path or die $!;
print {$default_out} encode_json($default_spec);
close $default_out or die $!;
is(system($^X, $driver, '--spec', $default_spec_path,
    '--step', 'plot_manhattan', '--list-steps'), 0,
    'default no-meta spec builds plot configuration');
my $default_runner = read_json("$dir/default_configs/auto_$default_spec->{artifact_stem}_runner.json");
is($default_runner->{TOP_HIT_MAF_THRESHOLD}, 0,
    'unknown GWAS MAF does not silently exclude every automatic top hit');
unlike($default_runner->{DISPLAY_GWAS}, qr/(?:^|,)\w*_?META(?:,|$)/,
    'default separate-input display does not contain a computed Meta track');
my $spec_path = "$dir/anemia.spec.json";
is(system($^X, $driver, '--gwas-dir', $dir, '--spec-out', $spec_path,
    '--generate-spec-only', '--no-compute-meta', '--no-exclude-strand-ambiguous',
    '--max-eaf-abs-diff', 'none', '--local-gtf-window-bp', '1e5'), 0,
    'anemia options generate a reusable spec');
my $spec = read_json($spec_path);
is($spec->{compute_meta}, 0, 'no-compute-meta survives spec generation');
is($spec->{exclude_strand_ambiguous}, 0, 'strand ambiguity opt-in survives spec generation');
is($spec->{max_eaf_abs_diff}, '', 'frequency threshold can be disabled explicitly');
is($spec->{local_gtf_window_bp}, '1e5', 'local window override survives spec generation');
is($spec->{exclude_local_gtf_tracks}, 'META', 'new local GTF specs omit meta track by default');
my %labels = map { $_->{label} => 1 } @{ $spec->{pairs} };
ok($labels{'EAS vs EUR+EAS'}, 'combined ancestry remains visible in pair label');
is($spec->{pairs}[0]{prefix}, 'EAS_EUR_PLUS_EAS', 'comparison prefix is safe for column names');
is($spec->{groups}[1]{tag}, 'GWAS_EUR_PLUS_EAS', 'combined ancestry gets a distinct group tag');
like(join(' ', @{ $spec->{auto_detect_summary}{warnings} || [] }), qr/overlapping samples|double count/i,
    'overlapping cohort comparison has a scientific warning');
my $sex_spec = { pairs => [{ group1_label => 'EUR_FEMALE', group2_label => 'EUR_MALE' }] };
is(cross_ancestry_filter_advisory($sex_spec), '',
    'female versus male within EUR is not mislabeled cross-ancestry');
$sex_spec->{pairs}[0]{group2_label} = 'EAS_MALE';
$sex_spec->{max_eaf_abs_diff} = 0.2;
$sex_spec->{top_hit_maf_threshold} = 0.01;
like(cross_ancestry_filter_advisory($sex_spec), qr/max_eaf_abs_diff=0\.2.*top_hit_maf_threshold=0\.01/,
    'explicit cutoffs in cross-ancestry comparison appear in the advisory');
$spec->{configs_dir} = "$dir/configs";
open my $spec_out, '>', $spec_path or die $!;
print {$spec_out} encode_json($spec);
close $spec_out or die $!;
is(system($^X, $driver, '--spec', $spec_path, '--step', 'diff_pairs', '--list-steps'), 0,
    'saved anemia spec prepares a differential run');
my $diff_config = read_json("$dir/configs/auto_$spec->{artifact_stem}_diff.json");
is($diff_config->{exclude_strand_ambiguous}, 0,
    'saved opt-in reaches the differential calculator');
is($diff_config->{max_eaf_abs_diff}, '',
    'disabled frequency threshold reaches the differential calculator');

my $catalog_dir = "$dir/catalog_inputs";
make_path($catalog_dir);
my $catalog_header = join("\t", qw(chromosome base_pair_location effect_allele other_allele
    beta standard_error effect_allele_frequency p_value variant_id hm_coordinate_conversion hm_code rsid));
for my $study (qw(GCST90018574 GCST90018794)) {
    my $table = "$catalog_header\n3\t108377251\tC\tG\t0.3\t0.2\t0.15\t0.15\tNA\tlo\t5\trsPalindromic\n";
    gzip \$table => "$catalog_dir/$study.h.tsv.gz" or die $GzipError;
}
my $catalog_spec_path = "$dir/catalog.spec.json";
is(system($^X, $driver, '--gwas-dir', $catalog_dir, '--spec-out', $catalog_spec_path,
    '--generate-spec-only'), 0,
    'original GWAS Catalog harmonised gzip headers generate a spec directly');
my $catalog_spec = read_json($catalog_spec_path);
is(scalar @{ $catalog_spec->{groups} }, 2, 'original accession files become separate GWAS groups');
like(join(' ', @{ $catalog_spec->{auto_detect_summary}{warnings} || [] }),
    qr/accession filenames do not identify ancestry/i,
    'accession-only grouping asks for study-aware ancestry labels');
like(join(' ', @{ $catalog_spec->{auto_detect_summary}{detected_files}[0]{resolved_columns} }),
    qr/POS=base_pair_location/, 'GWAS Catalog position column resolves without a rewritten header');
my $labeled_spec_path = "$dir/catalog.labeled.spec.json";
is(system($^X, $driver, '--gwas-dir', $catalog_dir, '--spec-out', $labeled_spec_path,
    '--generate-spec-only', '--group-label', 'GCST90018574=EAS',
    '--group-label', 'GCST90018794=EUR+EAS'), 0,
    'original accession files accept explicit ancestry labels');
my $labeled_spec = read_json($labeled_spec_path);
is($labeled_spec->{pairs}[0]{label}, 'EAS vs EUR+EAS',
    'original gzip input retains the correct cohort comparison label');
$labeled_spec->{configs_dir} = "$dir/catalog_configs";
open my $labeled_out, '>', $labeled_spec_path or die $!;
print {$labeled_out} encode_json($labeled_spec);
close $labeled_out or die $!;
is(system($^X, $driver, '--spec', $labeled_spec_path,
    '--step', 'plot_manhattan', '--list-steps'), 0,
    'original accession labels build a plotting configuration');
my $catalog_preset = read_json("$dir/catalog_configs/auto_$labeled_spec->{artifact_stem}_preset.json");
is($catalog_preset->{alias_map}{EUR_PLUS_EAS_P},
    'EAS_EUR_PLUS_EAS_GROUP2_P',
    'combined ancestry P column uses a valid plot identifier');
ok(!exists $catalog_preset->{alias_map}{'EUR+EAS_P'},
    'display punctuation does not enter plot data column names');
is(system($^X, $driver, '--spec', $labeled_spec_path, '--skip-plots',
    '--step', 'merge_raw', '--force'), 0,
    'original GWAS Catalog gzip files merge without header rewriting');
my $merged = '';
gunzip "$catalog_dir/$labeled_spec->{artifact_stem}_merged_long.tsv.gz" => \$merged
    or die $GunzipError;
like($merged, qr/rsPalindromic/, 'merged direct-import data retain the original rsID');
my $generated_diff = join("\t", qw(CHR BP A1 A2 SNP PAIR_TAG DIFF_BETA DIFF_SE DIFF_P)) . "\n"
    . "3\t108377251\tC\tG\trsPalindromic\tGWAS_EAS_vs_GWAS_EUR_PLUS_EAS\t0.1\t0.2\t0.6\n";
gzip \$generated_diff => "$catalog_dir/$labeled_spec->{artifact_stem}.tsv.gz" or die $GzipError;
my $rerun_spec_path = "$dir/catalog.rerun.spec.json";
is(system($^X, $driver, '--gwas-dir', $catalog_dir, '--spec-out', $rerun_spec_path,
    '--group-label', 'GCST90018574=EAS', '--group-label', 'GCST90018794=EUR+EAS',
    '--generate-spec-only'), 0,
    'directory rerun generates a spec with both original and generated files present');
is(read_json($rerun_spec_path)->{source_mode}, 'raw_pgc_vcf_sumstats',
    'directory rerun keeps the original cohort files as its source');

sub read_json {
    my ($path) = @_;
    open my $in, '<', $path or die $!;
    local $/;
    return decode_json(<$in>);
}

done_testing();
