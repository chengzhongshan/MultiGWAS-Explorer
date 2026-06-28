#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib $Bin;
use Getopt::Long qw(GetOptions);
use JSON::PP qw(decode_json);
use IO::Uncompress::Gunzip qw($GunzipError);
use File::Spec;
use File::Temp qw(tempfile);

my %opt = (
    top_hit_mode       => '',
    top_hit_focus_pvar => '',
    top_hit_signal_thrshd => '',
    top_hit_signal_thrshds => '',
    top_hit_signal_thrshd_fallback => '',
    top_hit_dist_bp    => '',
    maf_threshold      => '',
    gnomad_freq_file   => '',
    gnomad_pop_map     => '',
    gene_annotation_gtf => '',
    target_snps        => '',
    target_snp_genes   => '',
    max_hits           => 0,
    remove_x_chr       => 0,
);

GetOptions(
    'input=s'               => \$opt{input},
    'output=s'              => \$opt{output},
    'runner-config=s'       => \$opt{runner_config},
    'spec=s'                => \$opt{spec},
    'top-hit-mode=s'        => \$opt{top_hit_mode},
    'top-hit-focus-pvar=s'  => \$opt{top_hit_focus_pvar},
    'top-hit-signal-thrshd=s'  => \$opt{top_hit_signal_thrshd},
    'top-hit-signal-thrshds=s' => \$opt{top_hit_signal_thrshds},
    'top-hit-signal-thrshd-fallback=s' => \$opt{top_hit_signal_thrshd_fallback},
    'top-hit-dist-bp=s'     => \$opt{top_hit_dist_bp},
    'maf-threshold=f'       => \$opt{maf_threshold},
    'gnomad-freq-file=s'    => \$opt{gnomad_freq_file},
    'gnomad-pop-map=s'      => \$opt{gnomad_pop_map},
    'gene-annotation-gtf=s' => \$opt{gene_annotation_gtf},
    'target-snps=s'         => \$opt{target_snps},
    'target-snp-genes=s'    => \$opt{target_snp_genes},
    'max-hits=i'            => \$opt{max_hits},
    'remove-x-chr!'         => \$opt{remove_x_chr},
) or die usage();

die usage() unless length($opt{input} || '') && length($opt{output} || '');

my $runner = load_runner_config($opt{runner_config});
my $focus_pvar = first_nonempty(
    $opt{top_hit_focus_pvar},
    $runner->{TOP_HIT_FOCUS_PVAR},
    $runner->{MANHATTAN_P_VAR},
    'ALL_STD_P',
);
my $top_hit_mode = first_nonempty(
    $opt{top_hit_mode},
    $runner->{TOP_HIT_MODE},
    'differential',
);
my $thresholds = resolve_threshold_ladder(
    top_hit_mode        => $top_hit_mode,
    explicit_ladder     => $opt{top_hit_signal_thrshds},
    runner_ladder       => $runner->{TOP_HIT_SIGNAL_THRSHDS},
    primary_threshold   => first_nonempty(
        $opt{top_hit_signal_thrshd},
        $runner->{TOP_HIT_SIGNAL_THRSHD},
        '1e-6',
    ),
    fallback_threshold  => first_nonempty(
        $opt{top_hit_signal_thrshd_fallback},
        $runner->{TOP_HIT_SIGNAL_THRSHD_FALLBACK},
        '',
    ),
);
my $top_hit_dist_bp = first_nonempty(
    $opt{top_hit_dist_bp},
    $runner->{TOP_HIT_DIST_BP},
    '1e6',
);
my $maf_threshold = length($opt{maf_threshold} // '')
    ? $opt{maf_threshold}
    : first_nonempty($runner->{TOP_HIT_MAF_THRESHOLD}, 0.01);
my $gnomad_freq_file = first_nonempty(
    $opt{gnomad_freq_file},
    $runner->{TOP_HIT_GNOMAD_FREQ_FILE},
    '',
);
my $gnomad_pop_map = first_nonempty(
    $opt{gnomad_pop_map},
    $runner->{TOP_HIT_GNOMAD_POP_MAP},
    '',
);
my $gene_annotation_gtf = first_nonempty(
    $opt{gene_annotation_gtf},
    default_gene_annotation_gtf($runner->{REFERENCE_BUILD}),
    '',
);
my $target_snps = first_nonempty(
    $opt{target_snps},
    $runner->{TARGET_SNP_LIST},
    '',
);

my @hits;
if (length $target_snps) {
    my @targets = grep { length } map { trim($_) } split /,/, $target_snps;
    my %gene_map = parse_target_snp_gene_map($opt{target_snp_genes} || $runner->{TARGET_SNP_GENES});
    my $order = 0;
    @hits = map {
        $order++;
        +{
            SNP      => $_,
            hit_order => $order,
            gene     => $gene_map{ uc($_) } // '',
            snp_gene => length($gene_map{ uc($_) } // '') ? "$_:$gene_map{ uc($_) }" : '',
        }
    } @targets;
} elsif ($top_hit_mode =~ /^common_association$/i) {
    @hits = run_common_selector(
        input            => $opt{input},
        output           => $opt{output},
        runner_config    => $opt{runner_config},
        spec             => $opt{spec},
        thresholds       => $thresholds,
        top_hit_dist_bp  => $top_hit_dist_bp,
        max_hits         => $opt{max_hits},
        maf_threshold    => $maf_threshold,
        gnomad_freq_file => $gnomad_freq_file,
        gnomad_pop_map   => $gnomad_pop_map,
        remove_x_chr     => $opt{remove_x_chr},
    );
} elsif ($top_hit_mode =~ /^(?:common_and_differential|union_common_and_differential)$/i) {
    @hits = run_union_selector(
        input            => $opt{input},
        output           => $opt{output},
        runner_config    => $opt{runner_config},
        spec             => $opt{spec},
        focus_pvar       => $focus_pvar,
        thresholds       => $thresholds,
        top_hit_dist_bp  => $top_hit_dist_bp,
        max_hits         => $opt{max_hits},
        maf_threshold    => $maf_threshold,
        gnomad_freq_file => $gnomad_freq_file,
        gnomad_pop_map   => $gnomad_pop_map,
        remove_x_chr     => $opt{remove_x_chr},
    );
} else {
    @hits = run_differential_selector(
        input            => $opt{input},
        output           => $opt{output},
        focus_pvar       => $focus_pvar,
        thresholds       => $thresholds,
        top_hit_dist_bp  => $top_hit_dist_bp,
        max_hits         => $opt{max_hits},
        maf_threshold    => $maf_threshold,
        gnomad_freq_file => $gnomad_freq_file,
        gnomad_pop_map   => $gnomad_pop_map,
    );
}

die "No requested top hits were generated from $opt{input}\n" unless @hits;

my ($header, $rows) = build_export_rows(
    wide_data       => $opt{input},
    hits            => \@hits,
    focus_pvar      => $focus_pvar,
    target_snp_genes => ($opt{target_snp_genes} || $runner->{TARGET_SNP_GENES} || ''),
    gene_annotation_gtf => $gene_annotation_gtf,
);

write_csv(
    path   => $opt{output},
    header => $header,
    rows   => [ map { $rows->{ uc($_->{SNP} || '') } } @hits ],
);

print "Wrote requested top-hit CSV: $opt{output}\n";
print "Retained hits: " . scalar(@hits) . "\n";

sub run_common_selector {
    my (%args) = @_;
    my ($fh, $tmp) = tempfile('requested_top_hits_common_XXXXXX', SUFFIX => '.tsv', TMPDIR => 1, UNLINK => 1);
    close $fh;
    my @cmd = (
        $^X,
        File::Spec->catfile($Bin, 'verify_common_association_loci.pl'),
        '--input', $args{input},
        '--output', $tmp,
        '--top-p-thresholds', $args{thresholds},
        '--top-hit-dist-bp', $args{top_hit_dist_bp},
        '--max-loci', $args{max_hits},
        '--maf-threshold', $args{maf_threshold},
    );
    push @cmd, ('--runner-config', $args{runner_config})
        if defined $args{runner_config} && length $args{runner_config};
    push @cmd, ('--spec', $args{spec})
        if defined $args{spec} && length $args{spec};
    push @cmd, ('--gnomad-freq-file', $args{gnomad_freq_file})
        if defined $args{gnomad_freq_file} && length $args{gnomad_freq_file};
    push @cmd, ('--gnomad-pop-map', $args{gnomad_pop_map})
        if defined $args{gnomad_pop_map} && length $args{gnomad_pop_map};
    push @cmd, '--remove-x-chr' if $args{remove_x_chr};
    run_cmd(\@cmd, 'common-association selector');
    return read_hits_tsv($tmp);
}

sub run_differential_selector {
    my (%args) = @_;
    my ($fh, $tmp) = tempfile('requested_top_hits_diff_XXXXXX', SUFFIX => '.tsv', TMPDIR => 1, UNLINK => 1);
    close $fh;
    my $focus_prefix = ($args{focus_pvar} || '') =~ /^(.*?)_(?:STD_DIFF_P|STD_P|DIFF_P|GROUP1_P|GROUP2_P)$/
        ? $1
        : '';
    my @cmd = (
        $^X,
        File::Spec->catfile($Bin, 'gunplot', 'select_top_hits_from_wide.pl'),
        '--input', $args{input},
        '--output', $tmp,
        '--focus-pvar', $args{focus_pvar},
        '--thresholds', $args{thresholds},
        '--top-hit-dist-bp', $args{top_hit_dist_bp},
        '--max-hits', $args{max_hits},
        '--maf-threshold', $args{maf_threshold},
    );
    push @cmd, ('--focus-prefix', $focus_prefix) if length $focus_prefix;
    push @cmd, ('--gnomad-freq-file', $args{gnomad_freq_file})
        if defined $args{gnomad_freq_file} && length $args{gnomad_freq_file};
    push @cmd, ('--gnomad-pop-map', $args{gnomad_pop_map})
        if defined $args{gnomad_pop_map} && length $args{gnomad_pop_map};
    run_cmd(\@cmd, 'differential top-hit selector');
    return read_hits_tsv($tmp);
}

sub run_union_selector {
    my (%args) = @_;
    my @diff = run_differential_selector(%args);
    $_->{hit_class} ||= 'DIFFERENTIAL' for @diff;
    my @common = run_common_selector(%args);
    $_->{hit_class} ||= 'COMMON' for @common;
    return merge_hit_lists(
        max_hits => $args{max_hits},
        lists    => [ \@diff, \@common ],
    );
}

sub merge_hit_lists {
    my (%args) = @_;
    my %seen;
    my @merged;
    for my $list (@{ $args{lists} || [] }) {
        for my $hit (@{ $list || [] }) {
            next unless ref($hit) eq 'HASH';
            my $snp = uc($hit->{SNP} || '');
            next unless length $snp;
            next if $seen{$snp}++;
            push @merged, { %{$hit} };
        }
    }
    if (defined($args{max_hits}) && $args{max_hits} =~ /^\d+$/ && $args{max_hits} > 0 && @merged > $args{max_hits}) {
        @merged = @merged[0 .. $args{max_hits} - 1];
    }
    for my $i (0 .. $#merged) {
        $merged[$i]{hit_order} = $i + 1;
        $merged[$i]{panel_index} = $i + 1;
    }
    return @merged;
}

sub build_export_rows {
    my (%args) = @_;
    my %targets = map { uc($_->{SNP} || '') => 1 } grep { length($_->{SNP} || '') } @{ $args{hits} || [] };
    my %gene_override = parse_target_snp_gene_map($args{target_snp_genes});
    my %rows_by_snp;
    my (@wide_cols, %idx);

    my $fh = IO::Uncompress::Gunzip->new($args{wide_data})
        or die "Cannot open $args{wide_data}: $GunzipError\n";
    my $header = <$fh>;
    die "Wide data file is empty: $args{wide_data}\n" unless defined $header;
    chomp $header;
    $header =~ s/\r$//;
    @wide_cols = split /\t/, $header, -1;
    %idx = map { $wide_cols[$_] => $_ } 0 .. $#wide_cols;
    for my $required (qw(CHR BP SNP)) {
        die "Wide data file $args{wide_data} is missing required column $required\n"
            unless exists $idx{$required};
    }

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        my $snp = uc($f[$idx{SNP}] // '');
        next unless $targets{$snp};
        my %row;
        @row{@wide_cols} = @f;
        $rows_by_snp{$snp} = \%row;
        last if scalar(keys %rows_by_snp) == scalar(keys %targets);
    }
    close $fh;

    my %selected_chr = map {
        my $row = $rows_by_snp{ uc($_->{SNP} || '') } || {};
        my $chr = normalize_chr($row->{CHR} || $_->{CHR});
        length($chr) ? ($chr => 1) : ()
    } @{ $args{hits} || [] };
    my $genes_by_chr = {};
    if (defined $args{gene_annotation_gtf} && length($args{gene_annotation_gtf}) && -s $args{gene_annotation_gtf}) {
        $genes_by_chr = load_gtf_genes(
            file         => $args{gene_annotation_gtf},
            selected_chr => \%selected_chr,
        );
    }

    my @header = qw(
      hit_order panel_index CHR BP SNP EFFECT_ALLELE OTHER_ALLELE
      REFERENCE_ALLELE ALTERNATIVE_ALLELE gene snp_gene focus_signal hit_class
      selected_maf maf_source gwas_group1_maf gwas_group2_maf
      gwas_pair_maf_min gnomad_maf gnomad_pops maf_filter_decision
      maf_filter_reason
    );
    my %skip = map { $_ => 1 } qw(
      CHR BP SNP gene snp_gene A1 A2 REQUESTED_HIT_ORDER hit_order panel_index
      focus_signal EFFECT_ALLELE OTHER_ALLELE REFERENCE_ALLELE ALTERNATIVE_ALLELE
      selected_maf maf_source gwas_group1_maf gwas_group2_maf
      gwas_pair_maf_min gnomad_maf gnomad_pops maf_filter_decision
      maf_filter_reason
    );
    my @extra_cols = grep { !$skip{$_} } @wide_cols;
    push @header, @extra_cols;
    push @header, 'gene_source' unless grep { $_ eq 'gene_source' } @header;

    my %export_rows;
    my $fallback_order = 0;
    for my $hit (@{ $args{hits} || [] }) {
        $fallback_order++;
        my $snp = uc($hit->{SNP} || '');
        my %wide = %{ $rows_by_snp{$snp} || {} };
        next unless %wide;
        my %row = %wide;
        my $gene = length($gene_override{$snp} || '')
            ? $gene_override{$snp}
            : ($hit->{gene} || '');
        $gene = extract_gene_from_snp_gene($hit->{snp_gene}) if !length $gene;
        my $gene_source = defined $hit->{gene_source} && length($hit->{gene_source})
            ? $hit->{gene_source}
            : (defined $wide{gene_source} ? $wide{gene_source} : '');
        if (!length($gene) || $gene =~ /^(?:NA|N\/A|null)$/i) {
            my ($gtf_gene, $gtf_source) = resolve_gtf_gene_label(
                chr          => ($hit->{CHR} || $wide{CHR} || ''),
                bp           => ($hit->{BP} || $wide{BP} || ''),
                genes_by_chr => $genes_by_chr,
            );
            if (length($gtf_gene) && $gtf_gene !~ /^(?:NA|N\/A|null)$/i) {
                $gene = $gtf_gene;
                $gene_source = $gtf_source;
            }
        }
        my $snp_gene = length($gene_override{$snp} || '')
            ? ($hit->{SNP} || '') . ':' . $gene_override{$snp}
            : ($hit->{snp_gene} || '');
        $snp_gene = ($hit->{SNP} || '') . ':' . (length($gene) ? $gene : 'NA')
            unless length $snp_gene;
        my $focus_signal = (length($args{focus_pvar} || '') && exists $wide{ $args{focus_pvar} })
            ? $wide{ $args{focus_pvar} }
            : $hit->{focus_signal};
        my $hit_class = length($hit->{hit_class} || '')
            ? $hit->{hit_class}
            : normalize_hit_class($top_hit_mode);
        $row{hit_order} = defined $hit->{hit_order} && length($hit->{hit_order})
            ? $hit->{hit_order}
            : $fallback_order;
        $row{panel_index} = defined $hit->{panel_index} && length($hit->{panel_index})
            ? $hit->{panel_index}
            : $row{hit_order};
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
        $row{hit_class} = $hit_class;
        for my $maf_key (qw(selected_maf maf_source gwas_group1_maf gwas_group2_maf gwas_pair_maf_min gnomad_maf gnomad_pops maf_filter_decision maf_filter_reason)) {
            $row{$maf_key} = $hit->{$maf_key} if defined $hit->{$maf_key} && length $hit->{$maf_key};
        }
        $row{gene_source} = length($gene_source) ? $gene_source : (length($gene) ? 'NA' : 'NA');
        $export_rows{$snp} = \%row;
    }

    my @missing = grep { !exists $export_rows{$_} } sort keys %targets;
    if (@missing) {
        warn "WARNING: " . scalar(@missing) . " requested SNP(s) were not found in $args{wide_data}: "
            . join(', ', @missing) . "\n";
    }
    return (\@header, \%export_rows);
}

sub write_csv {
    my (%args) = @_;
    open my $fh, '>', $args{path} or die "Cannot write $args{path}: $!\n";
    print {$fh} join(',', map { csv_quote($_) } @{ $args{header} || [] }), "\n";
    for my $row (@{ $args{rows} || [] }) {
        next unless $row && ref($row) eq 'HASH';
        print {$fh} join(',', map { csv_quote($row->{$_}) } @{ $args{header} || [] }), "\n";
    }
    close $fh;
}

sub read_hits_tsv {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    my $header = <$fh>;
    die "Top-hit table is empty: $path\n" unless defined $header;
    chomp $header;
    $header =~ s/\r$//;
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
            panel_index  => pick_field(\%idx, \@f, qw(panel_index)),
            CHR          => pick_field(\%idx, \@f, qw(CHR chr)),
            BP           => pick_field(\%idx, \@f, qw(BP bp)),
            SNP          => pick_field(\%idx, \@f, qw(SNP snp rsid)),
            gene         => pick_field(\%idx, \@f, qw(gene nearest_gene genesymbol)),
            snp_gene     => pick_field(\%idx, \@f, qw(snp_gene SNP_GENE)),
            gene_source  => pick_field(\%idx, \@f, qw(gene_source)),
            focus_signal => pick_field(\%idx, \@f, qw(focus_signal common_assoc_p)),
            hit_class    => pick_field(\%idx, \@f, qw(hit_class HIT_CLASS hit_type HIT_TYPE)),
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
        (($a->{BP} || 0) <=> ($b->{BP} || 0))
    } @hits;
    return @hits;
}

sub parse_target_snp_gene_map {
    my ($text) = @_;
    my %map;
    return %map unless defined $text && length $text;
    for my $entry (split /,/, $text) {
        $entry = trim($entry);
        next unless length $entry;
        my ($snp, $gene) = split /:/, $entry, 2;
        $snp = trim($snp);
        $gene = trim($gene // '');
        next unless length $snp && length $gene;
        $map{ uc($snp) } = $gene;
    }
    return %map;
}

sub load_runner_config {
    my ($path) = @_;
    return {} unless defined $path && length $path;
    open my $fh, '<', $path or die "Cannot read runner config $path: $!\n";
    local $/;
    my $json = <$fh>;
    close $fh;
    my $cfg = decode_json($json);
    die "Runner config root must be a JSON object\n" unless ref($cfg) eq 'HASH';
    return $cfg;
}

sub run_cmd {
    my ($cmd, $label) = @_;
    my @cmd = @{$cmd || []};
    die "Internal error: empty command for $label\n" unless @cmd;
    print "[run] $label\n";
    system { $cmd[0] } @cmd;
    my $rc = $? >> 8;
    die "Command failed for $label (exit=$rc): @cmd\n" if $rc != 0;
}

sub pick_field {
    my ($idx, $fields, @candidates) = @_;
    for my $name (@candidates) {
        next unless exists $idx->{$name};
        my $value = $fields->[ $idx->{$name} ];
        return $value if defined $value && length $value;
    }
    return '';
}

sub normalize_hit_class {
    my ($mode) = @_;
    $mode = lc(trim($mode // ''));
    return 'COMMON' if $mode eq 'common_association';
    return 'COMMON' if $mode eq 'common_and_differential';
    return 'DIFFERENTIAL' if $mode eq 'differential';
    return 'SINGLE_GWAS' if $mode eq 'single_gwas';
    return 'CUSTOM_TARGET' if $mode eq 'targeted' || $mode eq 'custom_target';
    return 'CUSTOM_TARGET' unless length $mode;
    $mode =~ s/[^A-Za-z0-9]+/_/g;
    return uc($mode);
}

sub default_gene_annotation_gtf {
    my ($build) = @_;
    my $base = File::Spec->catdir($Bin, '..', 'cache', 'gtf');
    $build = lc(trim($build // 'hg38'));
    my @candidates;
    if ($build eq 'hg19' || $build eq 'grch37' || $build eq 'b37' || $build eq 'lift37') {
        @candidates = (
            File::Spec->catfile($base, 'gencode.v49lift37.annotation.gtf.gz'),
        );
    }
    elsif ($build eq 't2t' || $build eq 'hs1' || $build eq 'chm13' || $build eq 'chm13v2.0') {
        @candidates = (
            File::Spec->catfile($base, 'hs1.ncbiRefSeq.gtf.gz'),
        );
    }
    else {
        @candidates = (
            File::Spec->catfile($base, 'gencode.v49.annotation.gtf.gz'),
        );
    }
    for my $path (@candidates) {
        return $path if -s $path;
    }
    return '';
}

sub extract_gene_from_snp_gene {
    my ($value) = @_;
    return '' unless defined $value && length $value;
    my (undef, $gene) = split /:/, $value, 2;
    return trim($gene // '');
}

sub load_gtf_genes {
    my (%args) = @_;
    my $file = $args{file};
    my $selected_chr = $args{selected_chr} || {};
    my $fh = IO::Uncompress::Gunzip->new($file)
        or die "Cannot read GTF $file: $GunzipError\n";
    my %genes_by_chr;
    while (my $line = <$fh>) {
        next if $line =~ /^#/;
        chomp $line;
        $line =~ s/\r$//;
        my @f = split /\t/, $line, -1;
        next unless @f >= 9;
        next unless ($f[2] || '') eq 'gene';
        my $chr = normalize_chr($f[0]);
        next unless length $chr;
        next if %{$selected_chr} && !$selected_chr->{$chr};
        my $attrs = parse_gtf_attributes($f[8]);
        my $gene = $attrs->{gene_name} || $attrs->{gene_id} || '';
        next unless length $gene;
        push @{ $genes_by_chr{$chr} }, {
            gene  => $gene,
            start => 0 + $f[3],
            end   => 0 + $f[4],
            type  => ($attrs->{gene_type} || $attrs->{gene_biotype} || ''),
        };
    }
    close $fh;
    return \%genes_by_chr;
}

sub parse_gtf_attributes {
    my ($raw) = @_;
    my %attrs;
    while ($raw =~ /(\S+)\s+"([^"]*)"/g) {
        $attrs{$1} = $2;
    }
    return \%attrs;
}

sub resolve_gtf_gene_label {
    my (%args) = @_;
    my $chr = normalize_chr($args{chr});
    my $bp = pick_numeric($args{bp});
    return ('NA', 'NA') unless length($chr) && defined $bp;
    my $genes = $args{genes_by_chr}{$chr} || [];
    return ('NA', 'NA') unless @{$genes};

    my ($best_overlap_pc, $best_overlap_any, $best_nearest_pc, $best_nearest_any);
    for my $g (@{$genes}) {
        my $is_pc = (($g->{type} || '') eq 'protein_coding') ? 1 : 0;
        if ($bp >= $g->{start} && $bp <= $g->{end}) {
            if ($is_pc) {
                $best_overlap_pc = choose_better_gene($best_overlap_pc, $g, $bp);
            }
            else {
                $best_overlap_any = choose_better_gene($best_overlap_any, $g, $bp);
            }
        }
        else {
            if ($is_pc) {
                $best_nearest_pc = choose_better_gene($best_nearest_pc, $g, $bp);
            }
            else {
                $best_nearest_any = choose_better_gene($best_nearest_any, $g, $bp);
            }
        }
    }
    my $picked = $best_overlap_pc || $best_overlap_any || $best_nearest_pc || $best_nearest_any;
    return ('NA', 'NA') unless $picked;
    return ($picked->{gene}, 'GTF');
}

sub choose_better_gene {
    my ($cur, $cand, $bp) = @_;
    return $cand unless $cur;
    my $cur_dist = gene_distance($cur, $bp);
    my $cand_dist = gene_distance($cand, $bp);
    return $cand if $cand_dist < $cur_dist;
    return $cur if $cand_dist > $cur_dist;
    my $cur_span = ($cur->{end} || 0) - ($cur->{start} || 0);
    my $cand_span = ($cand->{end} || 0) - ($cand->{start} || 0);
    return $cand if $cand_span < $cur_span;
    return $cur if $cand_span > $cur_span;
    return (($cand->{gene} || '') cmp ($cur->{gene} || '')) < 0 ? $cand : $cur;
}

sub gene_distance {
    my ($g, $bp) = @_;
    return 0 if $bp >= $g->{start} && $bp <= $g->{end};
    return $g->{start} - $bp if $bp < $g->{start};
    return $bp - $g->{end};
}

sub pick_numeric {
    my ($value) = @_;
    return undef unless defined $value;
    $value = trim($value);
    return undef unless length $value;
    return undef unless $value =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
    return 0 + $value;
}

sub normalize_chr {
    my ($chr) = @_;
    return '' unless defined $chr;
    $chr = trim($chr);
    $chr =~ s/^chr//i;
    return 'X' if $chr =~ /^(?:23|X)$/i;
    return $chr if $chr =~ /^\d+$/;
    return '';
}

sub chr_order {
    my ($chr) = @_;
    return 9e9 unless defined $chr && length $chr;
    my $c = uc(trim($chr));
    $c =~ s/^CHR//;
    return 23 if $c eq 'X';
    return 24 if $c eq 'Y';
    return 25 if $c eq 'M' || $c eq 'MT';
    return $c =~ /^\d+$/ ? $c + 0 : 9e9;
}

sub csv_quote {
    my ($value) = @_;
    $value = '' unless defined $value;
    $value =~ s/\r\n?|\n/ /g;
    if ($value =~ /[",]/) {
        $value =~ s/"/""/g;
        return qq{"$value"};
    }
    return $value;
}

sub trim {
    my ($text) = @_;
    $text = '' unless defined $text;
    $text =~ s/^\s+//;
    $text =~ s/\s+$//;
    return $text;
}

sub first_nonempty {
    for my $value (@_) {
        next unless defined $value;
        return $value if ref($value) ? 1 : length("$value");
    }
    return '';
}

sub resolve_threshold_ladder {
    my (%args) = @_;
    my $explicit_ladder = trim($args{explicit_ladder});
    return $explicit_ladder if length $explicit_ladder;

    my $primary = trim($args{primary_threshold});
    $primary = '1e-6' unless length $primary;
    my $mode = lc(trim($args{top_hit_mode}));
    my $runner_ladder = trim($args{runner_ladder});
    my @runner_vals = unique_threshold_values(split /[,\s]+/, $runner_ladder);

    if ($mode ne 'differential') {
        return @runner_vals ? join(' ', @runner_vals) : $primary;
    }

    my $fallback = differential_threshold_fallback_value($primary, $args{fallback_threshold});
    my @vals = @runner_vals ? @runner_vals : ($primary);
    if (@vals <= 1 && length $fallback) {
        push @vals, $fallback unless grep { $_ eq $fallback } @vals;
    }
    return join(' ', unique_threshold_values(@vals));
}

sub differential_threshold_fallback_value {
    my ($primary, $fallback) = @_;
    $primary = trim($primary);
    $fallback = trim($fallback);
    return $fallback if length $fallback;
    return '1e-5' if normalized_numeric_text($primary) eq normalized_numeric_text('1e-6');
    return '';
}

sub unique_threshold_values {
    my @values = @_;
    my @uniq;
    for my $value (@values) {
        next unless defined $value && length trim($value);
        $value = trim($value);
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

sub usage {
    return <<"USAGE";
Usage:
  perl generate_requested_top_hits_csv.pl --input wide.tsv.gz --output top_hits.csv
      [--runner-config runner.json]
      [--top-hit-mode differential|common_association]
      [--top-hit-focus-pvar ALL_STD_P]
      [--top-hit-signal-thrshd 1e-6]
      [--top-hit-signal-thrshd-fallback 1e-5]
      [--target-snps rs1,rs2]
      [--target-snp-genes rs1:GENE1,rs2:GENE2]
      [--maf-threshold 0.01]
      [--gnomad-freq-file lookup.tsv.gz]
      [--gnomad-pop-map EUR:NFE,ASN:EAS]
USAGE
}
