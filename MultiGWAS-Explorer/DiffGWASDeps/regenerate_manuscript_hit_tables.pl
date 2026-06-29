#!/usr/bin/env perl
use strict;
use warnings;

use FindBin qw($Bin);
use lib $Bin;

use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw($GunzipError);
use JSON::PP qw(decode_json);
use Text::ParseWords qw(parse_line);
use TopHitMAF qw(
  numeric
  format_num
  derive_effect_af
  maf_from_effect_af
);

sub usage {
    return <<"USAGE";
Usage:
  perl DiffGWASDeps/regenerate_manuscript_hit_tables.pl [options]

Options:
  --config FILE                  Spec JSON used to infer the default wide table.
  --wide FILE                    Standardized wide GWAS table (.tsv.gz).
  --common-loci FILE             Validated common-locus TSV.
  --gtf FILE                     Gencode GTF.GZ used for fallback gene labels.
  --output-dir DIR               Output directory for manuscript tables.
  --representative-common-n N    Number of leading common loci for main Table 2.
  --diff-focus-pvar NAME         Differential focus P column. Default: ALL_STD_DIFF_P
  --diff-thresholds LIST         Differential threshold ladder. Default: 1e-6,1e-5
  --top-hit-dist-bp N            Distance-pruning span for differential loci. Default: 1e6
  --maf-threshold NUM            Minimum allowed differential MAF. Default: 0.01
USAGE
}

my %opt = (
    config                  => 'configs/spec_pgc_scz_sex_common_automation.json',
    common_loci             => 'tmp_common_verify_postfix_1e6.tsv',
    gtf                     => 'cache/gtf/gencode.v49lift37.annotation.gtf.gz',
    output_dir              => 'manuscript_assets/tables',
    representative_common_n => 15,
    diff_focus_pvar         => 'ALL_STD_DIFF_P',
    diff_thresholds         => '1e-6,1e-5',
    top_hit_dist_bp         => '1e6',
    maf_threshold           => 0.01,
);

GetOptions(
    'config=s'                  => \$opt{config},
    'wide=s'                    => \$opt{wide},
    'common-loci=s'             => \$opt{common_loci},
    'gtf=s'                     => \$opt{gtf},
    'output-dir=s'              => \$opt{output_dir},
    'representative-common-n=i' => \$opt{representative_common_n},
    'diff-focus-pvar=s'         => \$opt{diff_focus_pvar},
    'diff-thresholds=s'         => \$opt{diff_thresholds},
    'top-hit-dist-bp=s'         => \$opt{top_hit_dist_bp},
    'maf-threshold=f'           => \$opt{maf_threshold},
) or die usage();

$opt{wide} ||= infer_default_wide_from_config($opt{config});
die usage() unless $opt{wide};
die "Wide GWAS table not found: $opt{wide}\n" unless -s $opt{wide};
die "Common-locus TSV not found: $opt{common_loci}\n" unless -s $opt{common_loci};
die "GTF file not found: $opt{gtf}\n" unless -s $opt{gtf};
make_path($opt{output_dir}) unless -d $opt{output_dir};

my @diff_thresholds = grep { defined && length } map { trim($_) } split /[,\s]+/, ($opt{diff_thresholds} // '');
@diff_thresholds = ('1e-6', '1e-5') unless @diff_thresholds;
my $max_diff_threshold = max_num(map { 0 + $_ } @diff_thresholds);

my $common_loci = load_common_loci($opt{common_loci});
my %common_by_snp = map { uc($_->{SNP}) => $_ } @$common_loci;

my @reuse_label_files = grep { -s $_ } (
    catpath($opt{output_dir}, 'Table_2_representative_common_loci.csv'),
    catpath($opt{output_dir}, 'Table_S1_all_common_association_loci.csv'),
    catpath($opt{output_dir}, 'Table_S2_differential_loci.csv'),
);
my $gene_overrides = load_gene_overrides(@reuse_label_files);

my ($common_wide_rows, $diff_candidates) = scan_wide_for_hits(
    wide_file          => $opt{wide},
    common_by_snp      => \%common_by_snp,
    diff_focus_pvar    => $opt{diff_focus_pvar},
    max_diff_threshold => $max_diff_threshold,
    maf_threshold      => $opt{maf_threshold},
);

my ($diff_loci, $chosen_diff_threshold) = select_differential_loci(
    candidates       => $diff_candidates,
    thresholds       => \@diff_thresholds,
    top_hit_dist_bp  => $opt{top_hit_dist_bp},
);

my %selected_chr = map { normalize_chr($_->{CHR}) => 1 } (
    @$common_loci,
    @$diff_loci,
);
my $genes_by_chr = load_gtf_genes(
    file          => $opt{gtf},
    selected_chr  => \%selected_chr,
);

my @all_common_rows;
for my $loc (sort {
        ($a->{locus_rank} || 0) <=> ($b->{locus_rank} || 0)
    } @$common_loci) {
    my $wide = $common_wide_rows->{ uc($loc->{SNP}) } || {};
    my ($gene, $gene_source) = resolve_gene_label(
        snp            => $loc->{SNP},
        chr            => $loc->{CHR},
        bp             => $loc->{BP},
        overrides      => $gene_overrides,
        genes_by_chr   => $genes_by_chr,
    );

    push @all_common_rows, {
        hit_order         => $loc->{locus_rank},
        CHR               => $loc->{CHR},
        BP                => $loc->{BP},
        SNP               => $loc->{SNP},
        gene              => $gene,
        snp_gene          => build_snp_gene($loc->{SNP}, $gene),
        focus_signal      => $loc->{common_assoc_p},
        COMMON_ASSOC_P    => $loc->{common_assoc_p},
        ALL_STD_DIFF_P    => value_or_blank($wide->{ALL_STD_DIFF_P}),
        EUR_STD_DIFF_P    => value_or_blank($wide->{EUR_STD_DIFF_P}),
        ASN_STD_DIFF_P    => value_or_blank($wide->{ASN_STD_DIFF_P}),
        ALL_FEMALE_P      => value_or_blank($wide->{ALL_GROUP1_P}),
        ALL_FEMALE_BETA   => value_or_blank($wide->{ALL_GROUP1_BETA}),
        ALL_FEMALE_SE     => value_or_blank($wide->{ALL_GROUP1_SE}),
        ALL_MALE_P        => value_or_blank($wide->{ALL_GROUP2_P}),
        ALL_MALE_BETA     => value_or_blank($wide->{ALL_GROUP2_BETA}),
        ALL_MALE_SE       => value_or_blank($wide->{ALL_GROUP2_SE}),
        ALL_DIFF_P        => value_or_blank($wide->{ALL_DIFF_P}),
        ALL_DIFF_BETA     => value_or_blank($wide->{ALL_DIFF_BETA}),
        ALL_DIFF_SE       => value_or_blank($wide->{ALL_DIFF_SE}),
        EUR_FEMALE_P      => value_or_blank($wide->{EUR_GROUP1_P}),
        EUR_FEMALE_BETA   => value_or_blank($wide->{EUR_GROUP1_BETA}),
        EUR_FEMALE_SE     => value_or_blank($wide->{EUR_GROUP1_SE}),
        EUR_MALE_P        => value_or_blank($wide->{EUR_GROUP2_P}),
        EUR_MALE_BETA     => value_or_blank($wide->{EUR_GROUP2_BETA}),
        EUR_MALE_SE       => value_or_blank($wide->{EUR_GROUP2_SE}),
        EUR_DIFF_P        => value_or_blank($wide->{EUR_DIFF_P}),
        EUR_DIFF_BETA     => value_or_blank($wide->{EUR_DIFF_BETA}),
        EUR_DIFF_SE       => value_or_blank($wide->{EUR_DIFF_SE}),
        ASN_FEMALE_P      => value_or_blank($wide->{ASN_GROUP1_P}),
        ASN_FEMALE_BETA   => value_or_blank($wide->{ASN_GROUP1_BETA}),
        ASN_FEMALE_SE     => value_or_blank($wide->{ASN_GROUP1_SE}),
        ASN_MALE_P        => value_or_blank($wide->{ASN_GROUP2_P}),
        ASN_MALE_BETA     => value_or_blank($wide->{ASN_GROUP2_BETA}),
        ASN_MALE_SE       => value_or_blank($wide->{ASN_GROUP2_SE}),
        ASN_DIFF_P        => value_or_blank($wide->{ASN_DIFF_P}),
        ASN_DIFF_BETA     => value_or_blank($wide->{ASN_DIFF_BETA}),
        ASN_DIFF_SE       => value_or_blank($wide->{ASN_DIFF_SE}),
        selected_maf      => $loc->{selected_maf},
        maf_source        => $loc->{maf_source},
        gene_source       => $gene_source,
    };
}

my @representative_common_rows = @all_common_rows[ 0 .. min_num($#all_common_rows, $opt{representative_common_n} - 1) ];

my @diff_rows;
for my $loc (@$diff_loci) {
    my ($gene, $gene_source) = resolve_gene_label(
        snp            => $loc->{SNP},
        chr            => $loc->{CHR},
        bp             => $loc->{BP},
        overrides      => $gene_overrides,
        genes_by_chr   => $genes_by_chr,
    );
    push @diff_rows, {
        hit_order         => $loc->{hit_order},
        CHR               => $loc->{CHR},
        BP                => $loc->{BP},
        SNP               => $loc->{SNP},
        gene              => $gene,
        snp_gene          => build_snp_gene($loc->{SNP}, $gene),
        focus_signal      => $loc->{focus_signal},
        ALL_STD_DIFF_P    => $loc->{ALL_STD_DIFF_P},
        EUR_STD_DIFF_P    => value_or_blank($loc->{EUR_STD_DIFF_P}),
        ASN_STD_DIFF_P    => value_or_blank($loc->{ASN_STD_DIFF_P}),
        ALL_FEMALE_P      => value_or_blank($loc->{ALL_GROUP1_P}),
        ALL_MALE_P        => value_or_blank($loc->{ALL_GROUP2_P}),
        ALL_FEMALE_BETA   => value_or_blank($loc->{ALL_GROUP1_BETA}),
        ALL_FEMALE_SE     => value_or_blank($loc->{ALL_GROUP1_SE}),
        ALL_MALE_BETA     => value_or_blank($loc->{ALL_GROUP2_BETA}),
        ALL_MALE_SE       => value_or_blank($loc->{ALL_GROUP2_SE}),
        ALL_DIFF_P        => value_or_blank($loc->{ALL_DIFF_P}),
        ALL_DIFF_BETA     => value_or_blank($loc->{ALL_DIFF_BETA}),
        ALL_DIFF_SE       => value_or_blank($loc->{ALL_DIFF_SE}),
        EUR_FEMALE_P      => value_or_blank($loc->{EUR_GROUP1_P}),
        EUR_FEMALE_BETA   => value_or_blank($loc->{EUR_GROUP1_BETA}),
        EUR_FEMALE_SE     => value_or_blank($loc->{EUR_GROUP1_SE}),
        EUR_MALE_P        => value_or_blank($loc->{EUR_GROUP2_P}),
        EUR_MALE_BETA     => value_or_blank($loc->{EUR_GROUP2_BETA}),
        EUR_MALE_SE       => value_or_blank($loc->{EUR_GROUP2_SE}),
        EUR_DIFF_P        => value_or_blank($loc->{EUR_DIFF_P}),
        EUR_DIFF_BETA     => value_or_blank($loc->{EUR_DIFF_BETA}),
        EUR_DIFF_SE       => value_or_blank($loc->{EUR_DIFF_SE}),
        ASN_FEMALE_P      => value_or_blank($loc->{ASN_GROUP1_P}),
        ASN_FEMALE_BETA   => value_or_blank($loc->{ASN_GROUP1_BETA}),
        ASN_FEMALE_SE     => value_or_blank($loc->{ASN_GROUP1_SE}),
        ASN_MALE_P        => value_or_blank($loc->{ASN_GROUP2_P}),
        ASN_MALE_BETA     => value_or_blank($loc->{ASN_GROUP2_BETA}),
        ASN_MALE_SE       => value_or_blank($loc->{ASN_GROUP2_SE}),
        ASN_DIFF_P        => value_or_blank($loc->{ASN_DIFF_P}),
        ASN_DIFF_BETA     => value_or_blank($loc->{ASN_DIFF_BETA}),
        ASN_DIFF_SE       => value_or_blank($loc->{ASN_DIFF_SE}),
        selected_maf      => $loc->{selected_maf},
        maf_source        => $loc->{maf_source},
        gwas_group1_maf   => $loc->{gwas_group1_maf},
        gwas_group2_maf   => $loc->{gwas_group2_maf},
        gene_source       => $gene_source,
    };
}

write_csv(
    file   => catpath($opt{output_dir}, 'Table_S1_all_common_association_loci.csv'),
    header => [qw(
      hit_order CHR BP SNP gene snp_gene focus_signal COMMON_ASSOC_P
      ALL_STD_DIFF_P EUR_STD_DIFF_P ASN_STD_DIFF_P
      ALL_FEMALE_P ALL_FEMALE_BETA ALL_FEMALE_SE
      ALL_MALE_P ALL_MALE_BETA ALL_MALE_SE
      ALL_DIFF_P ALL_DIFF_BETA ALL_DIFF_SE
      EUR_FEMALE_P EUR_FEMALE_BETA EUR_FEMALE_SE
      EUR_MALE_P EUR_MALE_BETA EUR_MALE_SE
      EUR_DIFF_P EUR_DIFF_BETA EUR_DIFF_SE
      ASN_FEMALE_P ASN_FEMALE_BETA ASN_FEMALE_SE
      ASN_MALE_P ASN_MALE_BETA ASN_MALE_SE
      ASN_DIFF_P ASN_DIFF_BETA ASN_DIFF_SE
      selected_maf maf_source gene_source
    )],
    rows => \@all_common_rows,
);

write_csv(
    file   => catpath($opt{output_dir}, 'Table_S1_all_common_association_loci_full_strata.csv'),
    header => [qw(
      hit_order CHR BP SNP gene snp_gene focus_signal COMMON_ASSOC_P
      ALL_STD_DIFF_P EUR_STD_DIFF_P ASN_STD_DIFF_P
      ALL_FEMALE_P ALL_FEMALE_BETA ALL_FEMALE_SE
      ALL_MALE_P ALL_MALE_BETA ALL_MALE_SE
      ALL_DIFF_P ALL_DIFF_BETA ALL_DIFF_SE
      EUR_FEMALE_P EUR_FEMALE_BETA EUR_FEMALE_SE
      EUR_MALE_P EUR_MALE_BETA EUR_MALE_SE
      EUR_DIFF_P EUR_DIFF_BETA EUR_DIFF_SE
      ASN_FEMALE_P ASN_FEMALE_BETA ASN_FEMALE_SE
      ASN_MALE_P ASN_MALE_BETA ASN_MALE_SE
      ASN_DIFF_P ASN_DIFF_BETA ASN_DIFF_SE
      selected_maf maf_source gene_source
    )],
    rows => \@all_common_rows,
);

write_layout_csv(
    file      => catpath($opt{output_dir}, 'Table_S1_all_common_association_loci_xlsx_layout.csv'),
    row1      => [(
        'Locus summary', ('') x 4,
        'Standardized differential association', ('') x 2,
        'ALL_FEMALE', ('') x 2,
        'ALL_MALE', ('') x 2,
        'ALL_DIFF', ('') x 2,
        'EUR_FEMALE', ('') x 2,
        'EUR_MALE', ('') x 2,
        'EUR_DIFF', ('') x 2,
        'ASN_FEMALE', ('') x 2,
        'ASN_MALE', ('') x 2,
        'ASN_DIFF', ('') x 2,
    )],
    row2      => [qw(
        CHR BP SNP gene Smallest_ASSOC_P
        ALL_STD_DIFF_P EUR_STD_DIFF_P ASN_STD_DIFF_P
        BETA SE P
        BETA SE P
        BETA SE P
        BETA SE P
        BETA SE P
        BETA SE P
        BETA SE P
        BETA SE P
        BETA SE P
    )],
    src_cols   => [qw(
        CHR BP SNP gene COMMON_ASSOC_P
        ALL_STD_DIFF_P EUR_STD_DIFF_P ASN_STD_DIFF_P
        ALL_FEMALE_BETA ALL_FEMALE_SE ALL_FEMALE_P
        ALL_MALE_BETA ALL_MALE_SE ALL_MALE_P
        ALL_DIFF_BETA ALL_DIFF_SE ALL_DIFF_P
        EUR_FEMALE_BETA EUR_FEMALE_SE EUR_FEMALE_P
        EUR_MALE_BETA EUR_MALE_SE EUR_MALE_P
        EUR_DIFF_BETA EUR_DIFF_SE EUR_DIFF_P
        ASN_FEMALE_BETA ASN_FEMALE_SE ASN_FEMALE_P
        ASN_MALE_BETA ASN_MALE_SE ASN_MALE_P
        ASN_DIFF_BETA ASN_DIFF_SE ASN_DIFF_P
    )],
    rows      => \@all_common_rows,
);

write_csv(
    file   => catpath($opt{output_dir}, 'Table_2_representative_common_loci.csv'),
    header => [qw(
      hit_order CHR BP SNP gene focus_signal COMMON_ASSOC_P ALL_STD_DIFF_P
      ALL_FEMALE_BETA ALL_FEMALE_SE ALL_MALE_BETA ALL_MALE_SE ALL_DIFF_BETA ALL_DIFF_SE
      selected_maf maf_source
    )],
    rows => \@representative_common_rows,
);

write_csv(
    file   => catpath($opt{output_dir}, 'Table_S2_differential_loci.csv'),
    header => [qw(
      hit_order CHR BP SNP gene snp_gene focus_signal
      ALL_STD_DIFF_P EUR_STD_DIFF_P ASN_STD_DIFF_P
      ALL_FEMALE_P ALL_FEMALE_BETA ALL_FEMALE_SE
      ALL_MALE_P ALL_MALE_BETA ALL_MALE_SE
      ALL_DIFF_P ALL_DIFF_BETA ALL_DIFF_SE
      EUR_FEMALE_P EUR_FEMALE_BETA EUR_FEMALE_SE
      EUR_MALE_P EUR_MALE_BETA EUR_MALE_SE
      EUR_DIFF_P EUR_DIFF_BETA EUR_DIFF_SE
      ASN_FEMALE_P ASN_FEMALE_BETA ASN_FEMALE_SE
      ASN_MALE_P ASN_MALE_BETA ASN_MALE_SE
      ASN_DIFF_P ASN_DIFF_BETA ASN_DIFF_SE
      selected_maf maf_source gwas_group1_maf gwas_group2_maf gene_source
    )],
    rows => \@diff_rows,
);

write_csv(
    file   => catpath($opt{output_dir}, 'Table_S2_differential_loci_full_strata.csv'),
    header => [qw(
      hit_order CHR BP SNP gene snp_gene focus_signal
      ALL_STD_DIFF_P EUR_STD_DIFF_P ASN_STD_DIFF_P
      ALL_FEMALE_P ALL_FEMALE_BETA ALL_FEMALE_SE
      ALL_MALE_P ALL_MALE_BETA ALL_MALE_SE
      ALL_DIFF_P ALL_DIFF_BETA ALL_DIFF_SE
      EUR_FEMALE_P EUR_FEMALE_BETA EUR_FEMALE_SE
      EUR_MALE_P EUR_MALE_BETA EUR_MALE_SE
      EUR_DIFF_P EUR_DIFF_BETA EUR_DIFF_SE
      ASN_FEMALE_P ASN_FEMALE_BETA ASN_FEMALE_SE
      ASN_MALE_P ASN_MALE_BETA ASN_MALE_SE
      ASN_DIFF_P ASN_DIFF_BETA ASN_DIFF_SE
      selected_maf maf_source gwas_group1_maf gwas_group2_maf gene_source
    )],
    rows => \@diff_rows,
);

write_layout_csv(
    file      => catpath($opt{output_dir}, 'Table_S2_differential_loci_xlsx_layout.csv'),
    row1      => [(
        'Locus summary', ('') x 4,
        'Standardized differential association', ('') x 2,
        'ALL_FEMALE', ('') x 2,
        'ALL_MALE', ('') x 2,
        'ALL_DIFF', ('') x 2,
        'EUR_FEMALE', ('') x 2,
        'EUR_MALE', ('') x 2,
        'EUR_DIFF', ('') x 2,
        'ASN_FEMALE', ('') x 2,
        'ASN_MALE', ('') x 2,
        'ASN_DIFF', ('') x 2,
    )],
    row2      => [qw(
        CHR BP SNP gene Focus_DIFF_P
        ALL_STD_DIFF_P EUR_STD_DIFF_P ASN_STD_DIFF_P
        BETA SE P
        BETA SE P
        BETA SE P
        BETA SE P
        BETA SE P
        BETA SE P
        BETA SE P
        BETA SE P
        BETA SE P
    )],
    src_cols   => [qw(
        CHR BP SNP gene focus_signal
        ALL_STD_DIFF_P EUR_STD_DIFF_P ASN_STD_DIFF_P
        ALL_FEMALE_BETA ALL_FEMALE_SE ALL_FEMALE_P
        ALL_MALE_BETA ALL_MALE_SE ALL_MALE_P
        ALL_DIFF_BETA ALL_DIFF_SE ALL_DIFF_P
        EUR_FEMALE_BETA EUR_FEMALE_SE EUR_FEMALE_P
        EUR_MALE_BETA EUR_MALE_SE EUR_MALE_P
        EUR_DIFF_BETA EUR_DIFF_SE EUR_DIFF_P
        ASN_FEMALE_BETA ASN_FEMALE_SE ASN_FEMALE_P
        ASN_MALE_BETA ASN_MALE_SE ASN_MALE_P
        ASN_DIFF_BETA ASN_DIFF_SE ASN_DIFF_P
    )],
    rows      => \@diff_rows,
);

write_csv(
    file   => catpath($opt{output_dir}, 'Table_1_top_differential_locus.csv'),
    header => [qw(
      hit_order CHR BP SNP gene snp_gene focus_signal ALL_STD_DIFF_P ALL_FEMALE_P ALL_MALE_P
      ALL_FEMALE_BETA ALL_FEMALE_SE ALL_MALE_BETA ALL_MALE_SE ALL_DIFF_BETA ALL_DIFF_SE
      selected_maf maf_source gwas_group1_maf gwas_group2_maf gene_source
    )],
    rows => (@diff_rows ? [ $diff_rows[0] ] : []),
);

print "COMMON_LOCI\t" . scalar(@all_common_rows) . "\n";
print "REP_COMMON_LOCI\t" . scalar(@representative_common_rows) . "\n";
print "DIFF_CANDIDATES\t" . scalar(@$diff_candidates) . "\n";
print "DIFF_THRESHOLD\t" . ($chosen_diff_threshold // '') . "\n";
print "DIFF_LOCI\t" . scalar(@diff_rows) . "\n";
print "OUTPUT\t" . catpath($opt{output_dir}, 'Table_2_representative_common_loci.csv') . "\n";
print "OUTPUT\t" . catpath($opt{output_dir}, 'Table_S1_all_common_association_loci.csv') . "\n";
print "OUTPUT\t" . catpath($opt{output_dir}, 'Table_S1_all_common_association_loci_full_strata.csv') . "\n";
print "OUTPUT\t" . catpath($opt{output_dir}, 'Table_S1_all_common_association_loci_xlsx_layout.csv') . "\n";
print "OUTPUT\t" . catpath($opt{output_dir}, 'Table_S2_differential_loci.csv') . "\n";
print "OUTPUT\t" . catpath($opt{output_dir}, 'Table_S2_differential_loci_full_strata.csv') . "\n";
print "OUTPUT\t" . catpath($opt{output_dir}, 'Table_S2_differential_loci_xlsx_layout.csv') . "\n";
print "OUTPUT\t" . catpath($opt{output_dir}, 'Table_1_top_differential_locus.csv') . "\n";

sub infer_default_wide_from_config {
    my ($file) = @_;
    return '' unless defined $file && -s $file;
    local $/;
    open my $fh, '<', $file or die "Cannot read config $file: $!\n";
    my $json = <$fh>;
    close $fh;
    my $cfg = eval { decode_json($json) } || {};
    my $output_dir = native_path($cfg->{output_dir} || '');
    my $stem = $cfg->{artifact_stem} || '';
    return '' unless length($output_dir) && length($stem);
    my $sep = ($output_dir =~ m{[\\/]$}) ? '' : '/';
    return $output_dir . $sep . $stem . '.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz';
}

sub native_path {
    my ($path) = @_;
    return '' unless defined $path && length $path;
    if ($path =~ m{^/mnt/([A-Za-z])/(.*)$}) {
        my ($drive, $rest) = (uc($1), $2);
        $rest =~ s{/}{\\}g;
        return $drive . ':\\' . $rest;
    }
    return $path;
}

sub catpath {
    my ($dir, $file) = @_;
    return $dir =~ m{[\\/]$} ? $dir . $file : $dir . '/' . $file;
}

sub load_common_loci {
    my ($file) = @_;
    open my $fh, '<', $file or die "Cannot read common loci $file: $!\n";
    my $header = <$fh>;
    die "Common loci TSV is empty: $file\n" unless defined $header;
    chomp $header;
    $header =~ s/\r$//;
    my @cols = split /\t/, $header, -1;
    my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
    my @rows;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        push @rows, {
            locus_rank     => field_from_array(\@f, \%idx, 'locus_rank'),
            CHR            => field_from_array(\@f, \%idx, 'CHR'),
            BP             => field_from_array(\@f, \%idx, 'BP'),
            SNP            => field_from_array(\@f, \%idx, 'SNP'),
            common_assoc_p => field_from_array(\@f, \%idx, 'common_assoc_p'),
            selected_maf   => field_from_array(\@f, \%idx, 'selected_maf'),
            maf_source     => field_from_array(\@f, \%idx, 'maf_source'),
        };
    }
    close $fh;
    return \@rows;
}

sub load_gene_overrides {
    my @files = @_;
    my %map;
    for my $file (@files) {
        next unless -s $file;
        open my $fh, '<', $file or die "Cannot read gene override file $file: $!\n";
        my $header = <$fh>;
        next unless defined $header;
        chomp $header;
        $header =~ s/\r$//;
        my @cols = parse_csv($header);
        my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
        next unless exists $idx{SNP} && exists $idx{gene};
        while (my $line = <$fh>) {
            chomp $line;
            $line =~ s/\r$//;
            next unless length $line;
            my @f = parse_csv($line);
            my $snp = $f[ $idx{SNP} ] // '';
            my $gene = $f[ $idx{gene} ] // '';
            next unless length $snp;
            next unless length $gene;
            next if $gene =~ /^(?:NA|N\/A|null)$/i;
            my $key = uc($snp);
            next if exists $map{$key};
            $map{$key} = {
                gene        => $gene,
                gene_source => (exists $idx{gene_source} ? ($f[ $idx{gene_source} ] // '') : ''),
            };
        }
        close $fh;
    }
    return \%map;
}

sub scan_wide_for_hits {
    my (%args) = @_;
    my $wide_file = $args{wide_file};
    my $common_by_snp = $args{common_by_snp} || {};
    my $diff_focus_pvar = $args{diff_focus_pvar};
    my $max_diff_threshold = $args{max_diff_threshold};
    my $maf_threshold = defined $args{maf_threshold} ? $args{maf_threshold} : 0;

    my $fh = IO::Uncompress::Gunzip->new($wide_file)
        or die "Cannot read wide GWAS table $wide_file: $GunzipError\n";
    my $header = <$fh>;
    die "Wide GWAS table is empty: $wide_file\n" unless defined $header;
    chomp $header;
    $header =~ s/\r$//;
    my @cols = split /\t/, $header, -1;
    my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
    die "Missing required wide-table column: SNP\n" unless exists $idx{SNP};
    die "Missing required wide-table column: CHR\n" unless exists $idx{CHR};
    die "Missing required wide-table column: BP\n" unless exists $idx{BP};
    die "Missing required wide-table column: $diff_focus_pvar\n" unless exists $idx{$diff_focus_pvar};

    my %common_rows;
    my @diff_candidates;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        my $snp = field_from_array(\@f, \%idx, 'SNP');
        next unless length $snp;

        if (exists $common_by_snp->{ uc($snp) }) {
            $common_rows{ uc($snp) } = extract_wide_fields(\@f, \%idx);
        }

        my $diff_p = numeric(field_from_array(\@f, \%idx, $diff_focus_pvar));
        next unless defined $diff_p && $diff_p > 0 && $diff_p < $max_diff_threshold;
        my $diff_maf = annotate_diff_maf(
            raw           => \@f,
            idx           => \%idx,
            focus_prefix  => derive_focus_prefix($diff_focus_pvar),
            maf_threshold => $maf_threshold,
        );
        next unless $diff_maf->{pass};
        my $row = extract_wide_fields(\@f, \%idx);
        $row->{focus_signal} = format_num($diff_p);
        $row->{selected_maf} = format_num($diff_maf->{selected_maf});
        $row->{maf_source} = $diff_maf->{maf_source};
        $row->{gwas_group1_maf} = format_num($diff_maf->{gwas_group1_maf});
        $row->{gwas_group2_maf} = format_num($diff_maf->{gwas_group2_maf});
        push @diff_candidates, $row;
    }
    close $fh;

    return (\%common_rows, \@diff_candidates);
}

sub extract_wide_fields {
    my ($f, $idx) = @_;
    my %row;
    for my $name (qw(
        CHR BP SNP
        ALL_STD_DIFF_P EUR_STD_DIFF_P ASN_STD_DIFF_P
        ALL_GROUP1_P ALL_GROUP2_P ALL_DIFF_P
        ALL_GROUP1_BETA ALL_GROUP1_SE
        ALL_GROUP2_BETA ALL_GROUP2_SE
        ALL_DIFF_BETA ALL_DIFF_SE
        EUR_GROUP1_P EUR_GROUP2_P EUR_DIFF_P
        EUR_GROUP1_BETA EUR_GROUP1_SE
        EUR_GROUP2_BETA EUR_GROUP2_SE
        EUR_DIFF_BETA EUR_DIFF_SE
        ASN_GROUP1_P ASN_GROUP2_P ASN_DIFF_P
        ASN_GROUP1_BETA ASN_GROUP1_SE
        ASN_GROUP2_BETA ASN_GROUP2_SE
        ASN_DIFF_BETA ASN_DIFF_SE
    )) {
        $row{$name} = field_from_array($f, $idx, $name);
    }
    return \%row;
}

sub annotate_diff_maf {
    my (%args) = @_;
    my $g1_maf = compute_pair_group_maf($args{raw}, $args{idx}, $args{focus_prefix}, 1);
    my $g2_maf = compute_pair_group_maf($args{raw}, $args{idx}, $args{focus_prefix}, 2);
    my @vals = grep { defined $_ } ($g1_maf, $g2_maf);
    if (@vals) {
        @vals = sort { $a <=> $b } @vals;
        my $selected = $vals[0];
        return {
            selected_maf    => $selected,
            maf_source      => 'GWAS',
            gwas_group1_maf => $g1_maf,
            gwas_group2_maf => $g2_maf,
            pass            => ($args{maf_threshold} > 0 ? ($selected > $args{maf_threshold} ? 1 : 0) : 1),
        };
    }
    return {
        selected_maf    => undef,
        maf_source      => 'UNKNOWN',
        gwas_group1_maf => undef,
        gwas_group2_maf => undef,
        pass            => ($args{maf_threshold} > 0 ? 0 : 1),
    };
}

sub select_differential_loci {
    my (%args) = @_;
    my $candidates = $args{candidates} || [];
    my $thresholds = $args{thresholds} || [];
    my $top_hit_dist_bp = 0 + ($args{top_hit_dist_bp} || 0);
    my $half_window = $top_hit_dist_bp / 2;

    my $chosen_threshold;
    my @threshold_hits;
    for my $thr (@$thresholds) {
        my @cand = grep {
            my $p = numeric($_->{focus_signal});
            defined $p && $p > 0 && $p < (0 + $thr)
        } @$candidates;
        if (@cand) {
            $chosen_threshold = $thr;
            @threshold_hits = @cand;
            last;
        }
    }
    return ([], undef) unless @threshold_hits;

    @threshold_hits = sort {
        chr_order($a->{CHR}) <=> chr_order($b->{CHR})
            ||
        numeric($a->{focus_signal}) <=> numeric($b->{focus_signal})
            ||
        numeric($a->{BP}) <=> numeric($b->{BP})
    } @threshold_hits;

    my %selected_by_chr;
    my @selected;
    my $order = 0;
    for my $hit (@threshold_hits) {
        my $chr = normalize_chr($hit->{CHR});
        my $bp = numeric($hit->{BP});
        next unless length($chr) && defined $bp;
        my $keep = 1;
        for my $sel (@{ $selected_by_chr{$chr} || [] }) {
            my $sel_bp = numeric($sel->{BP});
            next unless defined $sel_bp;
            if ($bp >= $sel_bp - $half_window && $bp <= $sel_bp + $half_window) {
                $keep = 0;
                last;
            }
        }
        next unless $keep;
        $order++;
        $hit->{hit_order} = $order;
        push @selected, $hit;
        push @{ $selected_by_chr{$chr} }, $hit;
    }
    return (\@selected, $chosen_threshold);
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
        next if %$selected_chr && !$selected_chr->{$chr};
        my $attrs = parse_gtf_attributes($f[8]);
        my $gene = $attrs->{gene_name} || $attrs->{gene_id} || '';
        next unless length $gene;
        push @{ $genes_by_chr{$chr} }, {
            gene   => $gene,
            start  => 0 + $f[3],
            end    => 0 + $f[4],
            type   => ($attrs->{gene_type} || $attrs->{gene_biotype} || ''),
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

sub resolve_gene_label {
    my (%args) = @_;
    my $key = uc($args{snp} || '');
    if (length($key) && exists $args{overrides}{$key}) {
        my $hit = $args{overrides}{$key};
        return ($hit->{gene}, length($hit->{gene_source} || '') ? $hit->{gene_source} : 'MANUSCRIPT');
    }

    my $chr = normalize_chr($args{chr});
    my $bp = numeric($args{bp});
    return ('NA', 'NA') unless length($chr) && defined $bp;
    my $genes = $args{genes_by_chr}{$chr} || [];
    return ('NA', 'NA') unless @$genes;

    my ($best_overlap_pc, $best_overlap_any, $best_nearest_pc, $best_nearest_any);
    for my $g (@$genes) {
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

sub build_snp_gene {
    my ($snp, $gene) = @_;
    $gene = 'NA' unless defined $gene && length $gene;
    return ($snp || '') . ':' . $gene;
}

sub write_csv {
    my (%args) = @_;
    my $file = $args{file};
    my $header = $args{header} || [];
    my $rows = $args{rows} || [];
    open my $fh, '>', $file or die "Cannot write $file: $!\n";
    print {$fh} join(',', map { csv_quote($_) } @$header), "\n";
    for my $row (@$rows) {
        print {$fh} join(',', map { csv_quote(value_or_blank($row->{$_})) } @$header), "\n";
    }
    close $fh or die "Cannot close $file: $!\n";
}

sub write_layout_csv {
    my (%args) = @_;
    my $file = $args{file};
    my $row1 = $args{row1} || [];
    my $row2 = $args{row2} || [];
    my $src_cols = $args{src_cols} || [];
    my $rows = $args{rows} || [];

    die "layout csv row2/src_cols mismatch for $file\n"
        unless @$row2 == @$src_cols;
    die "layout csv row1/row2 mismatch for $file\n"
        unless @$row1 == @$row2;

    open my $fh, '>', $file or die "Cannot write $file: $!\n";
    print {$fh} join(',', map { csv_quote($_) } @$row1), "\n";
    print {$fh} join(',', map { csv_quote($_) } @$row2), "\n";
    for my $row (@$rows) {
        print {$fh} join(',', map { csv_quote(value_or_blank($row->{$_})) } @$src_cols), "\n";
    }
    close $fh or die "Cannot close $file: $!\n";
}

sub csv_quote {
    my ($v) = @_;
    $v = '' unless defined $v;
    $v =~ s/"/""/g;
    return qq("$v");
}

sub parse_csv {
    my ($line) = @_;
    my @row = parse_line(',', 0, $line);
    return map { defined $_ ? $_ : '' } @row;
}

sub field_from_array {
    my ($f, $idx, $name) = @_;
    return '' unless exists $idx->{$name};
    return trim($f->[ $idx->{$name} ] // '');
}

sub value_or_blank {
    my ($v) = @_;
    return '' unless defined $v;
    return $v;
}

sub compute_pair_group_maf {
    my ($raw, $idx, $prefix, $group_num) = @_;
    return undef unless defined $prefix && length $prefix;
    my $fa = numeric(field_from_array($raw, $idx, "${prefix}_GROUP${group_num}_FRQ_A"));
    my $fu = numeric(field_from_array($raw, $idx, "${prefix}_GROUP${group_num}_FRQ_U"));
    my $eaf = derive_effect_af(frq_a => $fa, frq_u => $fu);
    return maf_from_effect_af($eaf);
}

sub derive_focus_prefix {
    my ($focus_pvar) = @_;
    return '' unless defined $focus_pvar && length $focus_pvar;
    my $prefix = $focus_pvar;
    $prefix =~ s/_(?:STD_DIFF_P|STD_P|DIFF_P|GROUP1_P|GROUP2_P|P)$//;
    return $prefix;
}

sub normalize_chr {
    my ($chr) = @_;
    return '' unless defined $chr;
    $chr =~ s/^\s+|\s+$//g;
    $chr =~ s/^chr//i;
    return 'X' if $chr =~ /^(?:23|X)$/i;
    return $chr if $chr =~ /^\d+$/;
    return '';
}

sub chr_order {
    my ($chr) = @_;
    return 23 if defined $chr && uc($chr) eq 'X';
    return $chr if defined $chr && $chr =~ /^\d+$/;
    return 10_000;
}

sub min_num {
    my ($a, $b) = @_;
    return $a < $b ? $a : $b;
}

sub max_num {
    my (@vals) = @_;
    my $max = shift @vals;
    for my $v (@vals) {
        $max = $v if $v > $max;
    }
    return $max;
}

sub trim {
    my ($x) = @_;
    $x = '' unless defined $x;
    $x =~ s/^\s+//;
    $x =~ s/\s+$//;
    return $x;
}
