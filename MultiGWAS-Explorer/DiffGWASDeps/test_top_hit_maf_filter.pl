#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib $Bin;
use Getopt::Long qw(GetOptions);
use File::Spec;
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use IO::Compress::Gzip qw(gzip $GzipError);
use JSON::PP qw(decode_json encode_json);
use Text::ParseWords qw(parse_line);

my %opt = (
    workdir                => '',
    keep_workdir           => 0,
    no_real                => 0,
    real_common_thresholds => '5e-8 1e-6 1e-5',
);

GetOptions(
    'workdir=s'                => \$opt{workdir},
    'keep-workdir!'            => \$opt{keep_workdir},
    'no-real!'                 => \$opt{no_real},
    'real-wide=s'              => \$opt{real_wide},
    'real-runner-config=s'     => \$opt{real_runner_config},
    'real-common-thresholds=s' => \$opt{real_common_thresholds},
) or die usage();

my $workdir = length($opt{workdir})
    ? $opt{workdir}
    : tempdir('top_hit_maf_regress_XXXXXX', TMPDIR => 1, CLEANUP => ($opt{keep_workdir} ? 0 : 1));
make_path($workdir) unless -d $workdir;

print "Working directory: $workdir\n";

run_synthetic_differential_gwas_test($workdir);
run_synthetic_differential_threshold_fallback_test($workdir);
run_synthetic_differential_runner_ladder_fallback_test($workdir);
run_synthetic_differential_gnomad_test($workdir);
run_synthetic_common_test($workdir);

my ($real_wide, $real_runner_cfg) = resolve_real_inputs(\%opt);
if (!$opt{no_real} && defined $real_wide && defined $real_runner_cfg) {
    run_real_pgc_tests(
        workdir                => $workdir,
        real_wide              => $real_wide,
        real_runner_config     => $real_runner_cfg,
        real_common_thresholds => $opt{real_common_thresholds},
    );
}
else {
    print "Skipping real PGC schizophrenia validation because the required wide file or runner config was not available.\n";
}

print "All top-hit MAF safeguard regression checks passed.\n";
exit 0;

sub run_synthetic_differential_gwas_test {
    my ($dir) = @_;
    my $wide = File::Spec->catfile($dir, 'synthetic_diff_gwas_first.wide.tsv.gz');
    my $out  = File::Spec->catfile($dir, 'synthetic_diff_gwas_first.csv');
    write_tsv_gz(
        $wide,
        [qw(CHR BP SNP A1 A2 ALL_STD_P ALL_GROUP1_FRQ_A ALL_GROUP1_FRQ_U ALL_GROUP2_FRQ_A ALL_GROUP2_FRQ_U)],
        [
            [1, 100000, 'rsRareGWAS',   'A', 'G', '1e-9', '0.995', '0.995', '0.996', '0.996'],
            [1, 200000, 'rsCommonGWAS', 'C', 'T', '2e-9', '0.975', '0.975', '0.975', '0.975'],
        ],
    );
    run_cmd(
        [
            $^X,
            File::Spec->catfile($Bin, 'generate_requested_top_hits_csv.pl'),
            '--input', $wide,
            '--output', $out,
            '--top-hit-mode', 'differential',
            '--top-hit-focus-pvar', 'ALL_STD_P',
            '--top-hit-signal-thrshd', '1e-6',
            '--top-hit-signal-thrshds', '1e-6',
            '--top-hit-dist-bp', '1e8',
            '--maf-threshold', '0.01',
        ],
        'synthetic differential GWAS-first test',
    );
    my $rows = read_csv_rows($out);
    assert_eq(scalar(@$rows), 1, 'synthetic GWAS-first differential should retain exactly one top hit');
    my $row = $rows->[0];
    assert_eq($row->{SNP}, 'rsCommonGWAS', 'synthetic GWAS-first differential should keep rsCommonGWAS');
    assert_eq($row->{maf_source}, 'GWAS', 'synthetic GWAS-first differential should prefer GWAS MAF');
    assert_num_gt($row->{selected_maf}, 0.01, 'synthetic GWAS-first selected_maf should exceed threshold');
    assert_eq($row->{gwas_group1_maf}, '0.025', 'synthetic GWAS-first should export gwas_group1_maf');
    assert_eq($row->{gwas_group2_maf}, '0.025', 'synthetic GWAS-first should export gwas_group2_maf');
    print "PASS synthetic differential GWAS-first: SNP=$row->{SNP} selected_maf=$row->{selected_maf} source=$row->{maf_source}\n";
}

sub run_synthetic_differential_threshold_fallback_test {
    my ($dir) = @_;
    my $wide = File::Spec->catfile($dir, 'synthetic_diff_threshold_fallback.wide.tsv.gz');
    my $out  = File::Spec->catfile($dir, 'synthetic_diff_threshold_fallback.csv');
    write_tsv_gz(
        $wide,
        [qw(CHR BP SNP A1 A2 ALL_STD_P ALL_GROUP1_FRQ_A ALL_GROUP1_FRQ_U ALL_GROUP2_FRQ_A ALL_GROUP2_FRQ_U)],
        [
            [5, 100000, 'rsNeedsFallback', 'A', 'G', '4e-6', '0.97', '0.97', '0.975', '0.975'],
            [5, 200000, 'rsRareFallback',  'C', 'T', '3e-6', '0.995', '0.995', '0.996', '0.996'],
        ],
    );
    run_cmd(
        [
            $^X,
            File::Spec->catfile($Bin, 'generate_requested_top_hits_csv.pl'),
            '--input', $wide,
            '--output', $out,
            '--top-hit-mode', 'differential',
            '--top-hit-focus-pvar', 'ALL_STD_P',
            '--top-hit-signal-thrshd', '1e-6',
            '--top-hit-dist-bp', '1e8',
            '--maf-threshold', '0.01',
        ],
        'synthetic differential threshold fallback test',
    );
    my $rows = read_csv_rows($out);
    assert_eq(scalar(@$rows), 1, 'synthetic differential fallback should retain exactly one top hit');
    my $row = $rows->[0];
    assert_eq($row->{SNP}, 'rsNeedsFallback', 'synthetic differential fallback should keep rsNeedsFallback');
    assert_eq($row->{maf_source}, 'GWAS', 'synthetic differential fallback should still prefer GWAS MAF');
    assert_num_gt($row->{selected_maf}, 0.01, 'synthetic differential fallback selected_maf should exceed threshold');
    print "PASS synthetic differential threshold fallback: SNP=$row->{SNP} focus_signal=$row->{focus_signal} selected_maf=$row->{selected_maf}\n";
}

sub run_synthetic_differential_runner_ladder_fallback_test {
    my ($dir) = @_;
    my $wide = File::Spec->catfile($dir, 'synthetic_diff_runner_ladder_fallback.wide.tsv.gz');
    my $runner = File::Spec->catfile($dir, 'synthetic_diff_runner_ladder_fallback.json');
    my $out  = File::Spec->catfile($dir, 'synthetic_diff_runner_ladder_fallback.csv');
    write_tsv_gz(
        $wide,
        [qw(CHR BP SNP A1 A2 ALL_STD_P ALL_GROUP1_FRQ_A ALL_GROUP1_FRQ_U ALL_GROUP2_FRQ_A ALL_GROUP2_FRQ_U)],
        [
            [6, 100000, 'rsRunnerFallback', 'A', 'G', '6e-6', '0.97', '0.97', '0.975', '0.975'],
        ],
    );
    write_text(
        $runner,
        encode_json({
            TOP_HIT_MODE            => 'differential',
            TOP_HIT_FOCUS_PVAR      => 'ALL_STD_P',
            TOP_HIT_SIGNAL_THRSHD   => '1e-6',
            TOP_HIT_SIGNAL_THRSHDS  => '1e-6',
        })
    );
    run_cmd(
        [
            $^X,
            File::Spec->catfile($Bin, 'generate_requested_top_hits_csv.pl'),
            '--input', $wide,
            '--output', $out,
            '--runner-config', $runner,
            '--maf-threshold', '0.01',
        ],
        'synthetic differential runner-ladder fallback test',
    );
    my $rows = read_csv_rows($out);
    assert_eq(scalar(@$rows), 1, 'synthetic differential runner-ladder fallback should retain exactly one top hit');
    my $row = $rows->[0];
    assert_eq($row->{SNP}, 'rsRunnerFallback', 'synthetic differential runner-ladder fallback should keep rsRunnerFallback');
    assert_eq($row->{maf_source}, 'GWAS', 'synthetic differential runner-ladder fallback should prefer GWAS MAF');
    assert_num_gt($row->{selected_maf}, 0.01, 'synthetic differential runner-ladder fallback selected_maf should exceed threshold');
    print "PASS synthetic differential runner-ladder fallback: SNP=$row->{SNP} focus_signal=$row->{focus_signal} selected_maf=$row->{selected_maf}\n";
}

sub run_synthetic_differential_gnomad_test {
    my ($dir) = @_;
    my $wide = File::Spec->catfile($dir, 'synthetic_diff_gnomad_fallback.wide.tsv.gz');
    my $gnomad = File::Spec->catfile($dir, 'synthetic_gnomad_lookup.tsv');
    my $out  = File::Spec->catfile($dir, 'synthetic_diff_gnomad_fallback.csv');
    write_tsv_gz(
        $wide,
        [qw(CHR BP SNP A1 A2 EUR_STD_P)],
        [
            [4, 400000, 'rsCommonGNOMAD', 'T', 'C', '1e-9'],
        ],
    );
    write_text(
        $gnomad,
        join("\n",
            "SNP\tNFE_AF",
            "rsCommonGNOMAD\t0.02",
        ) . "\n"
    );
    run_cmd(
        [
            $^X,
            File::Spec->catfile($Bin, 'generate_requested_top_hits_csv.pl'),
            '--input', $wide,
            '--output', $out,
            '--top-hit-mode', 'differential',
            '--top-hit-focus-pvar', 'EUR_STD_P',
            '--top-hit-signal-thrshd', '1e-6',
            '--top-hit-signal-thrshds', '1e-6',
            '--top-hit-dist-bp', '1e8',
            '--maf-threshold', '0.01',
            '--gnomad-freq-file', $gnomad,
            '--gnomad-pop-map', 'EUR=NFE',
        ],
        'synthetic differential gnomAD-fallback test',
    );
    my $rows = read_csv_rows($out);
    assert_eq(scalar(@$rows), 1, 'synthetic gnomAD-fallback differential should retain exactly one top hit');
    my $row = $rows->[0];
    assert_eq($row->{SNP}, 'rsCommonGNOMAD', 'synthetic gnomAD-fallback differential should keep rsCommonGNOMAD');
    assert_eq($row->{maf_source}, 'GNOMAD', 'synthetic gnomAD-fallback differential should fall back to gnomAD');
    assert_num_gt($row->{selected_maf}, 0.01, 'synthetic gnomAD-fallback selected_maf should exceed threshold');
    assert_eq($row->{gnomad_maf}, '0.02', 'synthetic gnomAD-fallback should export gnomad_maf');
    print "PASS synthetic differential gnomAD-fallback: SNP=$row->{SNP} selected_maf=$row->{selected_maf} source=$row->{maf_source}\n";
}

sub run_synthetic_common_test {
    my ($dir) = @_;
    my $wide = File::Spec->catfile($dir, 'synthetic_common_association.wide.tsv.gz');
    my $runner = File::Spec->catfile($dir, 'synthetic_common_runner.json');
    my $out  = File::Spec->catfile($dir, 'synthetic_common_association.csv');
    write_tsv_gz(
        $wide,
        [qw(
            CHR BP A1 A2 SNP
            ALL_GROUP1_BETA ALL_GROUP2_BETA ALL_DIFF_BETA
            ALL_GROUP1_SE ALL_GROUP2_SE ALL_DIFF_SE
            ALL_GROUP1_P ALL_GROUP2_P ALL_DIFF_P ALL_STD_DIFF_Z ALL_STD_DIFF_P
            ALL_GROUP1_FRQ_A ALL_GROUP1_FRQ_U ALL_GROUP2_FRQ_A ALL_GROUP2_FRQ_U
            EUR_GROUP1_BETA EUR_GROUP2_BETA EUR_DIFF_BETA
            EUR_GROUP1_SE EUR_GROUP2_SE EUR_DIFF_SE
            EUR_GROUP1_P EUR_GROUP2_P EUR_DIFF_P EUR_STD_DIFF_Z EUR_STD_DIFF_P
            EUR_GROUP1_FRQ_A EUR_GROUP1_FRQ_U EUR_GROUP2_FRQ_A EUR_GROUP2_FRQ_U
        )],
        [
            [1,100000,'A','G','rsCommonCommon',0.5,0.02,0.48,0.05,0.10,0.11,'1e-8',0.4,'1e-6',4.1,'1e-4',0.97,0.98,0.97,0.98,0.3,0.01,0.29,0.09,0.11,0.14,0.03,0.8,0.2,1.1,0.2,0.97,0.98,0.97,0.98],
            [2,200000,'C','T','rsRareCommon',0.6,0.02,0.58,0.05,0.10,0.11,'5e-9',0.4,'1e-6',4.6,'1e-4',0.995,0.996,0.995,0.996,0.25,0.01,0.24,0.09,0.11,0.14,0.03,0.8,0.2,1.0,0.3,0.995,0.996,0.995,0.996],
        ],
    );
    write_text(
        $runner,
        encode_json({
            PAIR_DEFS => [
                { pair_tag => 'ALL_FEMALE_vs_MALE', group1 => 'ALL_FEMALE', group2 => 'ALL_MALE', prefix => 'ALL', label => 'All' },
                { pair_tag => 'EUR_FEMALE_vs_MALE', group1 => 'EUR_FEMALE', group2 => 'EUR_MALE', prefix => 'EUR', label => 'EUR' },
            ],
        })
    );
    run_cmd(
        [
            $^X,
            File::Spec->catfile($Bin, 'generate_requested_top_hits_csv.pl'),
            '--input', $wide,
            '--output', $out,
            '--top-hit-mode', 'common_association',
            '--runner-config', $runner,
            '--top-hit-signal-thrshds', '5e-8 1e-6 1e-5',
            '--top-hit-dist-bp', '1e8',
            '--maf-threshold', '0.01',
        ],
        'synthetic common-association MAF test',
    );
    my $rows = read_csv_rows($out);
    assert_eq(scalar(@$rows), 1, 'synthetic common-association test should retain exactly one top hit');
    my $row = $rows->[0];
    assert_eq($row->{SNP}, 'rsCommonCommon', 'synthetic common-association test should keep rsCommonCommon');
    assert_eq($row->{maf_source}, 'GWAS', 'synthetic common-association test should prefer GWAS MAF');
    assert_num_gt($row->{selected_maf}, 0.01, 'synthetic common-association selected_maf should exceed threshold');
    print "PASS synthetic common-association: SNP=$row->{SNP} selected_maf=$row->{selected_maf} source=$row->{maf_source}\n";
}

sub run_real_pgc_tests {
    my (%args) = @_;
    my $dir = $args{workdir};
    my $wide = $args{real_wide};
    my $runner = $args{real_runner_config};

    my $diff_out = File::Spec->catfile($dir, 'real_pgc_differential_top_hits.csv');
    run_cmd(
        [
            $^X,
            File::Spec->catfile($Bin, 'generate_requested_top_hits_csv.pl'),
            '--input', $wide,
            '--output', $diff_out,
            '--runner-config', $runner,
            '--maf-threshold', '0.01',
        ],
        'real PGC differential MAF validation',
    );
    my $diff_rows = read_csv_rows($diff_out);
    die "Real PGC differential validation returned no hits\n" unless @$diff_rows;
    for my $row (@$diff_rows) {
        assert_num_gt($row->{selected_maf}, 0.01, "real PGC differential hit $row->{SNP} should exceed MAF threshold");
        assert_true(length($row->{maf_source} || '') > 0, "real PGC differential hit $row->{SNP} should report maf_source");
    }
    my $lead_diff = $diff_rows->[0];
    print "PASS real PGC differential: lead SNP=$lead_diff->{SNP} selected_maf=$lead_diff->{selected_maf} source=$lead_diff->{maf_source} hits=" . scalar(@$diff_rows) . "\n";

    my $common_out = File::Spec->catfile($dir, 'real_pgc_common_top_hits.csv');
    run_cmd(
        [
            $^X,
            File::Spec->catfile($Bin, 'generate_requested_top_hits_csv.pl'),
            '--input', $wide,
            '--output', $common_out,
            '--runner-config', $runner,
            '--top-hit-mode', 'common_association',
            '--top-hit-signal-thrshds', $args{real_common_thresholds},
            '--maf-threshold', '0.01',
        ],
        'real PGC common-association MAF validation',
    );
    my $common_rows = read_csv_rows($common_out);
    die "Real PGC common-association validation returned no hits\n" unless @$common_rows;
    for my $row (@$common_rows) {
        assert_num_gt($row->{selected_maf}, 0.01, "real PGC common hit $row->{SNP} should exceed MAF threshold");
        assert_true(length($row->{maf_source} || '') > 0, "real PGC common hit $row->{SNP} should report maf_source");
    }
    my $lead_common = $common_rows->[0];
    print "PASS real PGC common-association: lead SNP=$lead_common->{SNP} selected_maf=$lead_common->{selected_maf} source=$lead_common->{maf_source} hits=" . scalar(@$common_rows) . "\n";
}

sub resolve_real_inputs {
    my ($opt) = @_;
    my $runner = $opt->{real_runner_config};
    my $wide = $opt->{real_wide};

    if (!defined $runner || !length $runner) {
        my $default_runner = File::Spec->catfile(File::Spec->updir(), 'configs', 'auto_PGC_SCZ_female_vs_male_diff_effects_runner.json');
        $default_runner = File::Spec->catfile($Bin, '..', 'configs', 'auto_PGC_SCZ_female_vs_male_diff_effects_runner.json');
        $runner = $default_runner if -f $default_runner;
    }
    return (undef, undef) unless defined $runner && -f $runner;

    if (!defined $wide || !length $wide) {
        my $cfg = load_json_file($runner);
        $wide = $cfg->{DATA_GZ} if ref($cfg) eq 'HASH' && defined $cfg->{DATA_GZ};
    }
    $wide = normalize_local_path($wide);
    return (-f $wide ? ($wide, $runner) : (undef, undef));
}

sub write_tsv_gz {
    my ($path, $header, $rows) = @_;
    my $content = join("\t", @$header) . "\n";
    for my $row (@$rows) {
        $content .= join("\t", @$row) . "\n";
    }
    gzip(\$content => $path) or die "Unable to write $path: $GzipError\n";
}

sub write_text {
    my ($path, $text) = @_;
    open my $fh, '>', $path or die "Cannot write $path: $!\n";
    print {$fh} $text;
    close $fh;
}

sub read_csv_rows {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    my $header = <$fh>;
    die "CSV is empty: $path\n" unless defined $header;
    chomp $header;
    my @header = parse_line(',', 1, $header);
    my @rows;
    while (my $line = <$fh>) {
        chomp $line;
        next unless length $line;
        my @vals = parse_line(',', 1, $line);
        my %row;
        @row{@header} = @vals;
        push @rows, \%row;
    }
    close $fh;
    return \@rows;
}

sub normalize_local_path {
    my ($path) = @_;
    return '' unless defined $path && length $path;
    if ($^O =~ /MSWin32/i) {
        $path =~ s{^/mnt/([A-Za-z])/$1:}{}i;
        $path =~ s{^/mnt/([A-Za-z])/(.*)$}{$1:/$2}i;
        $path =~ s{/}{\\}g if $path =~ m{^[A-Za-z]:/};
    }
    return $path;
}

sub load_json_file {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read JSON file $path: $!\n";
    local $/;
    my $json = <$fh>;
    close $fh;
    return decode_json($json);
}

sub run_cmd {
    my ($cmd, $label) = @_;
    print "[run] $label\n";
    print "       " . join(' ', map { shell_quote($_) } @$cmd) . "\n";
    my $rc = system(@$cmd);
    if ($rc != 0) {
        my $exit = $rc >> 8;
        die "Command failed for $label (exit=$exit)\n";
    }
}

sub shell_quote {
    my ($text) = @_;
    return "''" unless defined $text && length $text;
    if ($text =~ m{^[A-Za-z0-9_./:=-]+$}) {
        return $text;
    }
    $text =~ s/'/'"'"'/g;
    return "'$text'";
}

sub assert_true {
    my ($cond, $msg) = @_;
    die "ASSERTION FAILED: $msg\n" unless $cond;
}

sub assert_eq {
    my ($got, $want, $msg) = @_;
    $got  = '' unless defined $got;
    $want = '' unless defined $want;
    die "ASSERTION FAILED: $msg (got '$got', want '$want')\n" unless $got eq $want;
}

sub assert_num_gt {
    my ($got, $thr, $msg) = @_;
    die "ASSERTION FAILED: $msg (value missing)\n" unless defined $got && length $got;
    die "ASSERTION FAILED: $msg (got $got, expected > $thr)\n" unless 0 + $got > $thr;
}

sub usage {
    return <<"USAGE";
Usage: $0 [options]

Regression checks for the top-hit MAF safeguard.

Options:
  --workdir DIR                 Keep generated fixtures and outputs under DIR.
  --keep-workdir                Do not auto-delete the temporary workdir.
  --no-real                     Skip the real bundled PGC schizophrenia test.
  --real-wide PATH              Override the real bundled wide GWAS input.
  --real-runner-config PATH     Override the runner config used for the real test.
  --real-common-thresholds STR  Threshold ladder for real common-association validation.
                                Default: "5e-8 1e-6 1e-5"
USAGE
}
