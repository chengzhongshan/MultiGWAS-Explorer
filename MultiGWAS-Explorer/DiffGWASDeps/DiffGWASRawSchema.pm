package DiffGWASRawSchema;

use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw(basename);

our @EXPORT_OK = qw(
  normalize_header_name
  resolve_raw_header_aliases
  infer_raw_format_class
);

my %RAW_COLUMN_ALIASES = (
    CHROM   => [qw(CHROM CHR CHROMOSOME SEQNAME SEQID)],
    ID      => [qw(ID SNP RSID SNPID MARKER MARKERNAME VARIANTID VARIANT RS_NUMBER RS)],
    POS     => [qw(POS BP POSITION BASEPAIR BASEPAIRPOSITION BASE_PAIR_POSITION GENPOS PS)],
    A1      => [qw(A1 EA EFFECTALLELE EFFECT_ALLELE ALLELE1 TESTEDALLELE TESTED_ALLELE ALT INCALLELE CODEDALLELE)],
    A2      => [qw(A2 NEA NONEFFECTALLELE NON_EFFECT_ALLELE OTHERALLELE OTHER_ALLELE ALLELE2 REF DECALLELE NONCODEDALLELE ALLELE0)],
    FCAS    => [qw(FCAS FRQ_A CASEAF CASE_AF EAF_CASES AF_CASES FREQ_CASES CASE_FREQ CASE_FREQ_A)],
    FCON    => [qw(FCON FRQ_U CTRLAF CONTROLAF CONTROL_AF EAF_CONTROLS AF_CONTROLS FREQ_CONTROLS CONTROL_FREQ CTRL_FREQ)],
    IMPINFO => [qw(IMPINFO INFO INFO_SCORE IMPUTATIONINFO IMPUTATION_INFO RSQ MACHRSQ MACH_RSQ IMPUTEINFO)],
    BETA    => [qw(BETA EFFECT EFFECTSIZE EFFECT_SIZE B LOGOR LOG_ODDS ORBETA ESTIMATE BETA1 OR ODDSRATIO)],
    SE      => [qw(SE STDERR STANDARDERROR STANDARD_ERROR SEBETA BETA_SE)],
    PVAL    => [qw(PVAL PVALUE P_VALUE P PVALUELIN PVALMETA PV P_LRT PLRT P_WALD PWALD P_SCORE PSCORE)],
    NCAS    => [qw(NCAS N_CASES CASES CASECOUNT NCASE NCA TOT_CASES)],
    NCON    => [qw(NCON N_CONTROLS CONTROLS CONTROLCOUNT NCONTROL NCO TOT_CONTROLS)],
    NEFF    => [qw(NEFF N_EFFECTIVE EFFECTIVEN EFFECTIVE_N N_EFF)],
);

sub normalize_header_name {
    my ($text) = @_;
    $text = '' unless defined $text;
    $text = uc $text;
    $text =~ s/^\s+|\s+$//g;
    $text =~ s/[^A-Z0-9]+//g;
    return $text;
}

sub resolve_raw_header_aliases {
    my (%args) = @_;
    my $cols        = $args{cols}        || [];
    my $idx         = $args{idx}         || {};
    my $source_file = $args{source_file} // '';
    my $header      = $args{header}      // '';
    my $aliases     = merged_raw_column_aliases($args{alias_overrides});

    my %norm_to_actual;
    for my $col (@{$cols}) {
        my $norm = normalize_header_name($col);
        next unless length $norm;
        $norm_to_actual{$norm} ||= $col;
    }

    my %resolved;
    my @required = qw(CHROM ID POS A1 A2 BETA SE PVAL);
    my @optional = qw(FCAS FCON IMPINFO NCAS NCON NEFF);

    for my $canon (@required, @optional) {
        my $actual = '';
        for my $alias (@{ $aliases->{$canon} || [] }) {
            my $norm = normalize_header_name($alias);
            if (exists $norm_to_actual{$norm}) {
                $actual = $norm_to_actual{$norm};
                last;
            }
        }
        if (!length $actual) {
            $actual = find_special_raw_header_actual($canon, $cols);
        }
        if (length $actual) {
            $resolved{$canon} = $idx->{$actual};
        }
        elsif (grep { $_ eq $canon } @required) {
            die "Required column $canon not found in $source_file header.\nHeader was: $header\n";
        }
    }

    return %resolved;
}

sub infer_raw_format_class {
    my ($cols, $path) = @_;
    my %norm = map { normalize_header_name($_) => 1 } @{$cols || []};
    my $base = lc basename($path // '');
    return 'PGC_SUMSTATS_VCF' if $base =~ /\.vcf\.tsv\.gz$/i || ($norm{CHROM} && $norm{ID} && $norm{POS} && $norm{PVAL});
    return 'GEMMA_ASSOC' if $norm{CHR} && $norm{RS} && $norm{PS} && ($norm{PLRT} || $norm{PWALD} || $norm{PSCORE});
    return 'DANER_OR' if ($norm{OR} || $base =~ /^daner_/i);
    return 'TABULAR_GWAS';
}

sub find_special_raw_header_actual {
    my ($canon, $cols) = @_;
    for my $col (@{$cols}) {
        my $norm = normalize_header_name($col);
        return $col if $canon eq 'FCAS' && $norm =~ /^FRQA\d+$/;
        return $col if $canon eq 'FCON' && $norm =~ /^FRQU\d+$/;
    }
    return '';
}

sub merged_raw_column_aliases {
    my ($overrides) = @_;
    my %merged = map { $_ => [ @{ $RAW_COLUMN_ALIASES{$_} } ] } keys %RAW_COLUMN_ALIASES;
    return \%merged unless defined $overrides;
    die "raw_column_aliases must be a JSON object\n" unless ref($overrides) eq 'HASH';

    for my $raw_key (keys %{$overrides}) {
        my $canon = canonical_alias_key($raw_key);
        die "Unsupported raw_column_aliases key: $raw_key\n" unless exists $RAW_COLUMN_ALIASES{$canon};
        my $vals = $overrides->{$raw_key};
        my @extra = ref($vals) eq 'ARRAY' ? @{$vals} : ($vals);
        @extra = grep { defined $_ && length $_ } @extra;
        next unless @extra;
        my %seen;
        $merged{$canon} = [
            grep { !$seen{ normalize_header_name($_) }++ }
            (@extra, @{ $merged{$canon} || [] })
        ];
    }
    return \%merged;
}

sub canonical_alias_key {
    my ($key) = @_;
    my $norm = normalize_header_name($key);
    for my $canon (keys %RAW_COLUMN_ALIASES) {
        return $canon if normalize_header_name($canon) eq $norm;
    }
    return $norm;
}

1;
