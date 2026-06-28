#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib $Bin;
use Getopt::Long qw(GetOptions);
use JSON::PP qw(decode_json);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use File::Basename qw(basename);
use File::Spec;
use TopHitMAF qw(
  numeric
  derive_effect_af
  maf_from_effect_af
  parse_population_map
  infer_population_codes_for_text
  load_gnomad_lookup
  lookup_gnomad_maf
);

my $spec_file = '';
my $runner_config_file = '';
my $input_file = '';
my $output_file = '';
my $candidates_out = '';
my $top_p_threshold = '';
my $top_p_thresholds = '';
my $nominal_p_threshold = '0.05';
my $top_hit_dist_bp = '';
my $max_loci = 0;
my $direction_metric = 'beta';
my $remove_x_chr = 0;
my $maf_threshold = 0.01;
my $gnomad_freq_file = '';
my $gnomad_pop_map = '';
my $help = 0;

GetOptions(
    'spec=s'                => \$spec_file,
    'runner-config=s'       => \$runner_config_file,
    'input=s'               => \$input_file,
    'output=s'              => \$output_file,
    'candidates-out=s'      => \$candidates_out,
    'top-p-threshold=s'     => \$top_p_threshold,
    'top-p-thresholds=s'    => \$top_p_thresholds,
    'nominal-p-threshold=s' => \$nominal_p_threshold,
    'top-hit-dist-bp=s'     => \$top_hit_dist_bp,
    'max-loci=i'            => \$max_loci,
    'direction-metric=s'    => \$direction_metric,
    'remove-x-chr!'         => \$remove_x_chr,
    'maf-threshold=f'       => \$maf_threshold,
    'gnomad-freq-file=s'    => \$gnomad_freq_file,
    'gnomad-pop-map=s'      => \$gnomad_pop_map,
    'help!'                 => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

if (!$spec_file && !$runner_config_file && !$input_file) {
    die usage();
}

my $spec = {};
if ($spec_file) {
    $spec = load_json($spec_file);
}
my $runner_cfg = {};
if ($runner_config_file) {
    $runner_cfg = load_json($runner_config_file);
}

$input_file ||= infer_input_file($spec, $runner_cfg);
die "Unable to determine input wide GWAS file. Provide --input or a spec with artifact_stem/output_dir.\n"
    unless $input_file;

$input_file = normalize_local_path($input_file);
die "Input file not found: $input_file\n" unless -e $input_file;

my $artifact_stem = $spec->{artifact_stem} || $runner_cfg->{PROJECT_TAG} || guess_artifact_stem($input_file);
my $default_prefix = $artifact_stem . '.common_assoc_verify';
$output_file ||= File::Spec->catfile(File::Spec->curdir(), $default_prefix . '.tsv');
if (!$candidates_out) {
    if (defined $output_file && length $output_file) {
        my ($vol, $dir, $file) = File::Spec->splitpath($output_file);
        $file =~ s/\.tsv$//i;
        $candidates_out = File::Spec->catpath($vol, $dir, $file . '.candidates.tsv');
    }
    else {
        $candidates_out = File::Spec->catfile(File::Spec->curdir(), $default_prefix . '.candidates.tsv');
    }
}

my @pairs = infer_pairs($spec, $runner_cfg);
die "No pairs could be inferred from the spec or runner config. Provide a spec or runner config with pair/group tags.\n"
    unless @pairs;

my @thresholds = resolve_thresholds(
    top_p_threshold  => $top_p_threshold,
    top_p_thresholds => $top_p_thresholds,
);
$top_hit_dist_bp = length($top_hit_dist_bp)
    ? $top_hit_dist_bp
    : (defined $spec->{top_hit_dist_bp} ? $spec->{top_hit_dist_bp} : '1e6');

my $top_hit_dist_bp_num = 0 + $top_hit_dist_bp;
my $half_window_bp = $top_hit_dist_bp_num * 0.5;
my $nominal_p_threshold_num = 0 + $nominal_p_threshold;
my $pop_map = parse_population_map(
    length($gnomad_pop_map)
      ? $gnomad_pop_map
      : ($runner_cfg->{TOP_HIT_GNOMAD_POP_MAP} || $spec->{gnomad_population_map} || '')
);
my $resolved_gnomad_freq_file = length($gnomad_freq_file)
    ? $gnomad_freq_file
    : ($runner_cfg->{TOP_HIT_GNOMAD_FREQ_FILE} || $spec->{gnomad_freq_file} || '');
my $gnomad_lookup = load_gnomad_lookup(file => $resolved_gnomad_freq_file);

print STDERR "Reading wide GWAS table: $input_file\n";
print STDERR "Using top-hit threshold ladder: " . join(', ', @thresholds) . "\n";
print STDERR "Nominal partner threshold: $nominal_p_threshold_num\n";
print STDERR "Distance pruning window: total=$top_hit_dist_bp_num, half=$half_window_bp\n";
print STDERR "Direction metric: $direction_metric\n";
print STDERR "MAF threshold: $maf_threshold\n";

my ($fh, $close_fh) = open_table_handle($input_file);
my $header_line = <$fh>;
die "Input file is empty: $input_file\n" unless defined $header_line;
chomp $header_line;
$header_line =~ s/\r$//;
my @header = split /\t/, $header_line, -1;
my %idx = map { $header[$_] => $_ } 0 .. $#header;

for my $required (qw(CHR BP SNP)) {
    die "Missing required column '$required' in $input_file\n" unless exists $idx{$required};
}

for my $pair (@pairs) {
    $pair->{group1_pvar} = first_existing_column(\%idx, @{ $pair->{group1_pvar_candidates} });
    $pair->{group2_pvar} = first_existing_column(\%idx, @{ $pair->{group2_pvar_candidates} });
    $pair->{group1_betavar} = first_existing_column(\%idx, @{ $pair->{group1_betavar_candidates} });
    $pair->{group2_betavar} = first_existing_column(\%idx, @{ $pair->{group2_betavar_candidates} });
    $pair->{group1_zvar} = first_existing_column(\%idx, @{ $pair->{group1_zvar_candidates} });
    $pair->{group2_zvar} = first_existing_column(\%idx, @{ $pair->{group2_zvar_candidates} });
    $pair->{group1_frq_a_var} = first_existing_column(\%idx, @{ $pair->{group1_frq_a_candidates} || [] });
    $pair->{group1_frq_u_var} = first_existing_column(\%idx, @{ $pair->{group1_frq_u_candidates} || [] });
    $pair->{group2_frq_a_var} = first_existing_column(\%idx, @{ $pair->{group2_frq_a_candidates} || [] });
    $pair->{group2_frq_u_var} = first_existing_column(\%idx, @{ $pair->{group2_frq_u_candidates} || [] });

    for my $pvar ($pair->{group1_pvar}, $pair->{group2_pvar}) {
        die "Missing required P column '$pvar' in $input_file\n" unless exists $idx{$pvar};
    }
    for my $dirvar (grep { defined $_ && length $_ } ($pair->{group1_dirvar}, $pair->{group2_dirvar})) {
        die "Missing required direction column '$dirvar' in $input_file\n" unless exists $idx{$dirvar};
    }
}

my @groups = infer_groups_from_pairs(@pairs);
my %groups_by_name = map { $_->{name} => $_ } @groups;

my @candidates;
my %candidate_count_by_threshold;
my $maf_filtered_candidates = 0;
while (my $line = <$fh>) {
    chomp $line;
    $line =~ s/\r$//;
    next unless length $line;
    my @row = split /\t/, $line, -1;
    my $record = row_to_hash(\@header, \@row);
    next if $remove_x_chr && is_chr_x($record->{CHR});
    my ($common_p, @drivers) = compute_common_p($record, \@groups);
    next unless defined $common_p;

    my @qualifying_pairs = qualifying_common_groups(
        record      => $record,
        groups      => \@groups,
        common_p    => $common_p,
        nominal_thr => $nominal_p_threshold_num,
        direction_metric => $direction_metric,
    );
    next unless @qualifying_pairs;

    my $maf_info = common_candidate_maf_info(
        record         => $record,
        groups_by_name => \%groups_by_name,
        qualifying_pairs => \@qualifying_pairs,
        maf_threshold  => $maf_threshold,
        gnomad_lookup  => $gnomad_lookup,
        pop_map        => $pop_map,
    );
    if ($maf_threshold > 0 && !$maf_info->{pass}) {
        $maf_filtered_candidates++;
        next;
    }

    for my $thr (@thresholds) {
        next unless $common_p < $thr;
        $candidate_count_by_threshold{$thr}++;
    }

    push @candidates, {
        record           => $record,
        common_assoc_p   => $common_p,
        common_p_drivers => \@drivers,
        qualifying_pairs => \@qualifying_pairs,
        maf_info         => $maf_info,
    };
}
$close_fh->();

my $chosen_threshold;
for my $thr (@thresholds) {
    if (($candidate_count_by_threshold{$thr} || 0) > 0) {
        $chosen_threshold = $thr;
        last;
    }
}
$chosen_threshold = $thresholds[-1] unless defined $chosen_threshold;

my @threshold_hits = grep { $_->{common_assoc_p} < $chosen_threshold } @candidates;

print STDERR "Raw qualifying candidates across all thresholds: " . scalar(@candidates) . "\n";
for my $thr (@thresholds) {
    print STDERR "  candidates with COMMON_ASSOC_P < $thr : " . ($candidate_count_by_threshold{$thr} || 0) . "\n";
}
print STDERR "Chosen threshold for verification: $chosen_threshold\n";
print STDERR "Candidates surviving chosen threshold before distance pruning: " . scalar(@threshold_hits) . "\n";
print STDERR "Candidates filtered by MAF: $maf_filtered_candidates\n";

@threshold_hits = sort {
    $a->{record}{CHR} <=> $b->{record}{CHR}
        ||
    $a->{common_assoc_p} <=> $b->{common_assoc_p}
        ||
    $a->{record}{BP} <=> $b->{record}{BP}
} @threshold_hits;

my @selected;
my %selected_by_chr;
for my $hit (@threshold_hits) {
    my $chr = 0 + ($hit->{record}{CHR} // 0);
    my $bp  = 0 + ($hit->{record}{BP}  // 0);
    my $keep = 1;
    for my $sel (@{ $selected_by_chr{$chr} || [] }) {
        if ($bp >= $sel->{dis_st} && $bp <= $sel->{dis_end}) {
            $keep = 0;
            last;
        }
    }
    next unless $keep;

    $hit->{dis_st} = $bp - $half_window_bp;
    $hit->{dis_end} = $bp + $half_window_bp;
    push @selected, $hit;
    push @{ $selected_by_chr{$chr} }, $hit;
    last if $max_loci && @selected >= $max_loci;
}

print STDERR "Selected loci after distance pruning: " . scalar(@selected) . "\n";

write_output_table(
    path    => $candidates_out,
    header  => [qw(candidate_rank CHR BP SNP A1 A2 common_assoc_p common_p_drivers qualifying_pairs selected_maf maf_source gwas_pair_maf_min gnomad_maf gnomad_pops maf_filter_decision maf_filter_reason key)],
    rows    => [
        map {
            my $r = $_->{record};
            my $m = $_->{maf_info} || {};
            [
                undef,
                (map { value_or_blank($r->{$_}) } qw(CHR BP SNP A1 A2)),
                format_num($_->{common_assoc_p}),
                join(';', @{ $_->{common_p_drivers} }),
                join(';', @{ $_->{qualifying_pairs} }),
                format_num($m->{selected_maf}),
                value_or_blank($m->{maf_source}),
                format_num($m->{gwas_pair_maf_min}),
                format_num($m->{gnomad_maf}),
                value_or_blank($m->{gnomad_pops}),
                value_or_blank($m->{decision}),
                value_or_blank($m->{reason}),
                make_key($r),
            ]
        } @threshold_hits
    ],
    rank_col => 0,
);

write_output_table(
    path    => $output_file,
    header  => [qw(locus_rank CHR BP SNP A1 A2 common_assoc_p common_p_drivers qualifying_pairs selected_maf maf_source gwas_pair_maf_min gnomad_maf gnomad_pops maf_filter_decision maf_filter_reason key dis_st dis_end)],
    rows    => [
        map {
            my $r = $_->{record};
            my $m = $_->{maf_info} || {};
            [
                undef,
                (map { value_or_blank($r->{$_}) } qw(CHR BP SNP A1 A2)),
                format_num($_->{common_assoc_p}),
                join(';', @{ $_->{common_p_drivers} }),
                join(';', @{ $_->{qualifying_pairs} }),
                format_num($m->{selected_maf}),
                value_or_blank($m->{maf_source}),
                format_num($m->{gwas_pair_maf_min}),
                format_num($m->{gnomad_maf}),
                value_or_blank($m->{gnomad_pops}),
                value_or_blank($m->{decision}),
                value_or_blank($m->{reason}),
                make_key($r),
                $_->{dis_st},
                $_->{dis_end},
            ]
        } @selected
    ],
    rank_col => 0,
);

print STDERR "Wrote candidate table: $candidates_out\n";
print STDERR "Wrote pruned loci table: $output_file\n";

exit 0;

sub infer_input_file {
    my ($spec, $runner_cfg) = @_;
    if ($runner_cfg && ref($runner_cfg) eq 'HASH' && $runner_cfg->{DATA_GZ}) {
        return $runner_cfg->{DATA_GZ};
    }
    return '' unless $spec && ref($spec) eq 'HASH';

    if ($spec->{wide_input}) {
        return $spec->{wide_input};
    }

    my $stem = $spec->{artifact_stem} || return '';
    my $dir = $spec->{output_dir} || $spec->{input_dir} || return '';
    return join_path($dir, $stem . '.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz');
}

sub is_chr_x {
    my ($chr) = @_;
    return 0 unless defined $chr;
    $chr =~ s/^\s+|\s+$//g;
    $chr =~ s/^chr//i;
    return $chr =~ /^(?:X|23)$/i ? 1 : 0;
}

sub infer_pairs {
    my ($spec, $runner_cfg) = @_;
    my @pairs;
    my $pair_source =
        ($spec && ref($spec) eq 'HASH' && ref($spec->{pairs}) eq 'ARRAY') ? $spec->{pairs}
      : ($runner_cfg && ref($runner_cfg) eq 'HASH' && ref($runner_cfg->{PAIR_DEFS}) eq 'ARRAY') ? $runner_cfg->{PAIR_DEFS}
      : [];
    return @pairs unless ref($pair_source) eq 'ARRAY' && @{$pair_source};

    for my $pair (@{$pair_source}) {
        next unless ref($pair) eq 'HASH';
        my $g1 = $pair->{group1} || next;
        my $g2 = $pair->{group2} || next;
        my $prefix = $pair->{prefix} || '';
        my $g1_label = safe_name($pair->{group1_label} || $g1);
        my $g2_label = safe_name($pair->{group2_label} || $g2);
        push @pairs, {
            pair_tag     => ($pair->{pair_tag} || join('_vs_', $g1, $g2)),
            prefix       => $prefix,
            group1       => $g1,
            group2       => $g2,
            group1_pvar_candidates => [ $g1_label . '_P', ($prefix ? $prefix . '_GROUP1_P' : ()) ],
            group2_pvar_candidates => [ $g2_label . '_P', ($prefix ? $prefix . '_GROUP2_P' : ()) ],
            group1_betavar_candidates => [ $g1_label . '_BETA', ($prefix ? $prefix . '_GROUP1_BETA' : ()) ],
            group2_betavar_candidates => [ $g2_label . '_BETA', ($prefix ? $prefix . '_GROUP2_BETA' : ()) ],
            group1_zvar_candidates => [ $g1_label . '_Z', ($prefix ? $prefix . '_GROUP1_Z' : ()) ],
            group2_zvar_candidates => [ $g2_label . '_Z', ($prefix ? $prefix . '_GROUP2_Z' : ()) ],
            group1_frq_a_candidates => [ ($prefix ? $prefix . '_GROUP1_FRQ_A' : ()) ],
            group1_frq_u_candidates => [ ($prefix ? $prefix . '_GROUP1_FRQ_U' : ()) ],
            group2_frq_a_candidates => [ ($prefix ? $prefix . '_GROUP2_FRQ_A' : ()) ],
            group2_frq_u_candidates => [ ($prefix ? $prefix . '_GROUP2_FRQ_U' : ()) ],
        };
    }
    return @pairs;
}

sub resolve_thresholds {
    my (%args) = @_;
    if (length($args{top_p_threshold} || '')) {
        return (0 + $args{top_p_threshold});
    }
    if (length($args{top_p_thresholds} || '')) {
        return map { 0 + $_ } grep { length } split /[\s,]+/, $args{top_p_thresholds};
    }
    return (5e-8, 1e-6, 1e-5);
}

sub compute_common_p {
    my ($record, $groups) = @_;
    my %seen;
    my @pvars;
    for my $group (@$groups) {
        push @pvars, $group->{pvar};
    }
    @pvars = grep { !$seen{$_}++ } @pvars;

    my $min_p;
    for my $pvar (@pvars) {
        my $v = normalize_num($record->{$pvar});
        next unless defined $v && $v > 0;
        $min_p = $v if !defined($min_p) || $v < $min_p;
    }
    return unless defined $min_p;

    my @drivers = grep {
        my $v = normalize_num($record->{$_});
        defined($v) && $v > 0 && $v == $min_p
    } @pvars;

    return ($min_p, @drivers);
}

sub qualifying_common_groups {
    my (%args) = @_;
    my $record      = $args{record};
    my $groups      = $args{groups};
    my $common_p    = $args{common_p};
    my $nominal_thr = $args{nominal_thr};
    my $direction_metric = lc($args{direction_metric} || 'beta');

    my @qualified;
    for my $driver (@$groups) {
        my $p1 = normalize_num($record->{ $driver->{pvar} });
        next unless defined $p1 && $p1 > 0 && $p1 == $common_p;
        my $d1 = group_direction_value($record, $driver, $direction_metric);
        next unless defined $d1 && $d1 != 0;
        for my $partner (@$groups) {
            next if $partner->{name} eq $driver->{name};
            my $p2 = normalize_num($record->{ $partner->{pvar} });
            next unless defined $p2 && $p2 > 0 && $p2 <= $nominal_thr;
            my $d2 = group_direction_value($record, $partner, $direction_metric);
            next unless same_direction($d1, $d2);
            push @qualified, "$driver->{name}->$partner->{name}";
        }
    }
    return @qualified;
}

sub open_table_handle {
    my ($path) = @_;
    my $fh;
    if ($path =~ /\.gz$/i) {
        $fh = IO::Uncompress::Gunzip->new($path)
            or die "Unable to open gzip file $path: $GunzipError\n";
    }
    else {
        open($fh, '<', $path) or die "Unable to open $path: $!\n";
    }
    return ($fh, sub { close $fh; });
}

sub row_to_hash {
    my ($header, $row) = @_;
    my %h;
    @h{@$header} = @$row;
    return \%h;
}

sub write_output_table {
    my (%args) = @_;
    my $path = $args{path};
    my @header = @{ $args{header} };
    my @rows = @{ $args{rows} };
    my $rank_col = $args{rank_col};

    open(my $out, '>', $path) or die "Unable to write $path: $!\n";
    print {$out} join("\t", @header) . "\n";
    my $rank = 0;
    for my $row (@rows) {
        $rank++;
        $row->[$rank_col] = $rank if defined $rank_col;
        print {$out} join("\t", map { defined $_ ? $_ : '' } @$row) . "\n";
    }
    close $out;
}

sub load_json {
    my ($path) = @_;
    open(my $fh, '<', $path) or die "Unable to open $path: $!\n";
    local $/;
    my $json = <$fh>;
    close $fh;
    return decode_json($json);
}

sub normalize_local_path {
    my ($path) = @_;
    return '' unless defined $path;
    if ($path =~ m{^/mnt/([a-zA-Z])/(.*)$}) {
        my ($drive, $rest) = ($1, $2);
        $rest =~ s{/}{\\}g;
        return uc($drive) . ":\\" . $rest;
    }
    return $path;
}

sub join_path {
    my ($dir, $leaf) = @_;
    $dir =~ s{[\\/]+$}{};
    return $dir . '/' . $leaf;
}

sub guess_artifact_stem {
    my ($path) = @_;
    $path =~ s{\\}{/}g if defined $path;
    my $base = basename($path);
    $base =~ s/\.stdized\.wide_beta_se_p_p_lt_0p05\.final\.tsv\.gz$//;
    $base =~ s/\.tsv\.gz$//;
    return $base;
}

sub normalize_num {
    my ($v) = @_;
    return numeric($v);
}

sub format_num {
    my ($v) = @_;
    return '' unless defined $v;
    return sprintf('%.12g', $v);
}

sub group_direction_value {
    my ($record, $group, $metric) = @_;

    if ($metric eq 'z') {
        my $z = normalize_num($record->{ $group->{zvar} });
        return $z if defined $z;
    }

    my $beta = normalize_num($record->{ $group->{betavar} });
    return $beta if defined $beta;

    my $z = normalize_num($record->{ $group->{zvar} });
    return $z if defined $z;

    return undef;
}

sub infer_groups_from_pairs {
    my (@pairs) = @_;
    my %groups;
    for my $pair (@pairs) {
        $groups{$pair->{group1}} ||= {
            name    => $pair->{group1},
            pvar    => $pair->{group1_pvar},
            betavar => $pair->{group1_betavar},
            zvar    => $pair->{group1_zvar},
            frq_a_var => $pair->{group1_frq_a_var},
            frq_u_var => $pair->{group1_frq_u_var},
        };
        $groups{$pair->{group2}} ||= {
            name    => $pair->{group2},
            pvar    => $pair->{group2_pvar},
            betavar => $pair->{group2_betavar},
            zvar    => $pair->{group2_zvar},
            frq_a_var => $pair->{group2_frq_a_var},
            frq_u_var => $pair->{group2_frq_u_var},
        };
    }
    return sort { $a->{name} cmp $b->{name} } values %groups;
}

sub common_candidate_maf_info {
    my (%args) = @_;
    my $record = $args{record} || {};
    my $groups_by_name = $args{groups_by_name} || {};
    my $qualifying_pairs = $args{qualifying_pairs} || [];
    my $maf_threshold = defined $args{maf_threshold} ? $args{maf_threshold} : 0;

    my @passing;
    for my $qp (@{$qualifying_pairs}) {
        my ($driver, $partner) = split /->/, $qp, 2;
        next unless defined $driver && defined $partner;
        my $driver_group = $groups_by_name->{$driver} || next;
        my $partner_group = $groups_by_name->{$partner} || next;

        my $driver_maf = group_maf_from_record($record, $driver_group);
        my $partner_maf = group_maf_from_record($record, $partner_group);
        my @gwas_mafs = grep { defined $_ } ($driver_maf, $partner_maf);
        my $pair_maf_min;
        my %info = (
            pair => $qp,
        );
        if (@gwas_mafs) {
            @gwas_mafs = sort { $a <=> $b } @gwas_mafs;
            $pair_maf_min = $gwas_mafs[0];
            %info = (
                %info,
                selected_maf      => $pair_maf_min,
                maf_source        => 'GWAS',
                gwas_pair_maf_min => $pair_maf_min,
            );
        }
        else {
            my @pop_codes = infer_population_codes_for_text($driver, $args{pop_map});
            push @pop_codes, infer_population_codes_for_text($partner, $args{pop_map});
            my $gnomad = lookup_gnomad_maf(
                lookup    => $args{gnomad_lookup},
                record    => $record,
                pop_codes => \@pop_codes,
            );
            next unless $gnomad && defined $gnomad->{maf};
            my @pairs = map { $_ . '=' . format_num($gnomad->{pop_mafs}{$_}) } sort keys %{ $gnomad->{pop_mafs} || {} };
            %info = (
                %info,
                selected_maf => $gnomad->{maf},
                maf_source   => 'GNOMAD',
                gnomad_maf   => $gnomad->{maf},
                gnomad_pops  => join('|', @pairs),
            );
        }
        next if $maf_threshold > 0 && (!defined $info{selected_maf} || $info{selected_maf} <= $maf_threshold);
        push @passing, \%info;
    }

    if (@passing) {
        @passing = sort {
            ($b->{selected_maf} // -1) <=> ($a->{selected_maf} // -1)
              ||
            ($a->{pair} cmp $b->{pair})
        } @passing;
        my $best = $passing[0];
        return {
            selected_maf      => $best->{selected_maf},
            maf_source        => $best->{maf_source},
            gwas_pair_maf_min => $best->{gwas_pair_maf_min},
            gnomad_maf        => $best->{gnomad_maf},
            gnomad_pops       => $best->{gnomad_pops},
            decision          => 'PASS',
            reason            => sprintf('At least one qualifying pair exceeded the MAF threshold: %s', $best->{pair}),
            pass              => 1,
        };
    }

    return {
        maf_source => 'UNKNOWN',
        decision   => ($maf_threshold > 0 ? 'FILTERED' : 'PASS'),
        reason     => 'No qualifying driver/partner pair had usable MAF support above the threshold',
        pass       => ($maf_threshold > 0 ? 0 : 1),
    };
}

sub group_maf_from_record {
    my ($record, $group) = @_;
    return undef unless $record && $group;
    my $eaf = derive_effect_af(
        frq_a => $record->{ $group->{frq_a_var} },
        frq_u => $record->{ $group->{frq_u_var} },
    );
    return maf_from_effect_af($eaf);
}

sub same_direction {
    my ($a, $b) = @_;
    return 0 unless defined $a && defined $b;
    return 0 if $a == 0 || $b == 0;
    return (($a > 0 && $b > 0) || ($a < 0 && $b < 0)) ? 1 : 0;
}

sub value_or_blank {
    my ($v) = @_;
    return defined $v ? $v : '';
}

sub make_key {
    my ($r) = @_;
    if (defined $r->{Key} && length $r->{Key}) {
        return $r->{Key};
    }
    return join(':', map { value_or_blank($r->{$_}) } qw(CHR BP));
}

sub safe_name {
    my ($text) = @_;
    $text //= '';
    $text =~ s/[^A-Za-z0-9]+/_/g;
    $text =~ s/^_+|_+$//g;
    return uc($text);
}

sub first_existing_column {
    my ($idx, @candidates) = @_;
    for my $candidate (@candidates) {
        next unless defined $candidate && length $candidate;
        return $candidate if exists $idx->{$candidate};
    }
    return $candidates[0] // '';
}

sub usage {
    return <<'USAGE';
Usage:
  perl DiffGWASDeps/verify_common_association_loci.pl --spec spec.json
  perl DiffGWASDeps/verify_common_association_loci.pl --spec spec.json --top-p-threshold 5e-8
  perl DiffGWASDeps/verify_common_association_loci.pl --input wide.tsv.gz --spec spec.json --top-hit-dist-bp 5000000

Description:
  Independently verify "common association" loci from the wide differential GWAS
  table used by auto_prepare_and_run_diff_gwas.pl.

  The script:
    1. reads the wide GWAS table
    2. computes COMMON_ASSOC_P as the minimum single-GWAS P among all GWAS groups
    3. keeps rows where the best single-GWAS signal also has at least one other
       GWAS with nominal association in the same effect direction
    4. applies chromosome-wise distance pruning that matches
       get_top_signal_within_dist.sas:
         dis_st = BP - top_hit_dist_bp * 0.5
         dis_end = BP + top_hit_dist_bp * 0.5

Options:
  --spec FILE                Automation spec JSON with pairs/group tags.
  --input FILE               Override the inferred wide GWAS input file.
  --output FILE              Output TSV for distance-pruned loci.
  --candidates-out FILE      Output TSV for all candidates before pruning.
  --top-p-threshold NUM      Use a single common-association threshold.
  --top-p-thresholds LIST    Space/comma-separated threshold ladder.
  --nominal-p-threshold NUM  Partner nominal P threshold. Default: 0.05
  --top-hit-dist-bp NUM      Total exclusion window. Default: spec or 1e6
  --max-loci N               Stop after N selected loci.
  --direction-metric STR     beta or z. Default: beta
  --remove-x-chr             Exclude chromosome X / 23 from selected loci.
  --help                     Show this help.
USAGE
}
