#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $input_dir =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs';
my $output =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_sex_stratified_merged_long.tsv.gz';
my $manifest =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_sex_stratified_merged_long.manifest.tsv';
my $limit = 0;

GetOptions(
    'input-dir=s' => \$input_dir,
    'output=s'    => \$output,
    'manifest=s'  => \$manifest,
    'limit=i'     => \$limit,
) or die "Usage: $0 [--input-dir DIR] [--output OUT.tsv.gz] [--manifest MANIFEST.tsv] [--limit N]\n";

my @groups = (
    {
        tag   => 'SCZ_W3_ASN_FEMALE',
        files => [
            'daner_PGC_SCZ_w3_14_0618a_asn_female.gz',
            'daner_scz_w3_HRC_chrX_asn_fem_run2.gz',
        ],
    },
    {
        tag   => 'SCZ_W3_ASN_MALE',
        files => [
            'daner_PGC_SCZ_w3_14_0618a_asn_male.gz',
            'daner_scz_w3_HRC_chrX_asn_mal_run2.gz',
        ],
    },
    {
        tag   => 'SCZ_W3_EUR_FEMALE',
        files => [
            'daner_PGC_SCZ_w3_75_0618a_eur_female.gz',
            'daner_scz_w3_HRC_chrX_eur_fem_deduped_0518e.gz',
        ],
    },
    {
        tag   => 'SCZ_W3_EUR_MALE',
        files => [
            'daner_PGC_SCZ_w3_75_0618a_eur_male.gz',
            'daner_scz_w3_HRC_chrX_eur_mal_deduped_0518e.gz',
        ],
    },
    {
        tag   => 'SCZ_W3_ALL_FEMALE_AUTOSOME',
        files => ['daner_PGC_SCZ_w3_81_0618a_all_female.gz'],
    },
    {
        tag   => 'SCZ_W3_ALL_MALE_AUTOSOME',
        files => ['daner_PGC_SCZ_w3_81_0618a_all_male.gz'],
    },
    {
        tag   => 'SCZ_W3_UKBBDEDUPE_AUTOSOME',
        files => ['daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.gz'],
    },
);

my @out_cols = qw(
  CHR BP A1 A2 SNP GWAS_TAG SOURCE_FILE CHR_ORIGINAL IS_CHRX
  FRQ_A FRQ_U INFO OR SE P ngt Direction HetISqt HetDf HetPVa Nca Nco Neff Neff_half
);

open my $out_fh, '-|', 'true' or die "Internal pipe check failed: $!";
close $out_fh;

open my $out, '|-', "gzip -c > '$output'"
  or die "Cannot write gzip output to $output: $!\n";
print {$out} join("\t", @out_cols), "\n";

open my $man, '>', $manifest
  or die "Cannot write manifest to $manifest: $!\n";
print {$man} join("\t", qw(GWAS_TAG SOURCE_FILE ROWS_WRITTEN HEADER STATUS)), "\n";

for my $group (@groups) {
    for my $file (@{ $group->{files} }) {
        my $path = "$input_dir/$file";
        unless (-e $path) {
            print {$man} join("\t", $group->{tag}, $file, 0, '', 'MISSING'), "\n";
            warn "Missing expected file: $path\n";
            next;
        }

        my $rows = process_file($path, $file, $group->{tag}, $out, $limit);
        print {$man} join("\t", $group->{tag}, $file, $rows->{rows}, $rows->{header}, $rows->{status}), "\n";
        warn "Finished $file as $group->{tag}: $rows->{rows} rows\n";
    }
}

close $out or die "Failed closing gzip output $output: $!\n";
close $man or die "Failed closing manifest $manifest: $!\n";

warn "Merged long-format table: $output\n";
warn "Manifest: $manifest\n";

sub process_file {
    my ($path, $source_file, $tag, $out, $limit) = @_;

    open my $fh, '-|', "zcat '$path'"
      or die "Cannot read $path with zcat: $!\n";

    my $header = <$fh>;
    unless (defined $header) {
        close $fh;
        return { rows => 0, header => '', status => 'EMPTY_OR_UNREADABLE' };
    }
    chomp $header;
    $header =~ s/\r$//;

    my @cols = split /\t/, $header, -1;
    @cols = split /\s+/, $header if @cols == 1;
    my %idx;
    for my $i (0 .. $#cols) {
        $idx{$cols[$i]} = $i;
    }

    for my $required (qw(CHR SNP BP A1 A2 INFO OR SE P ngt Direction HetISqt HetDf HetPVa Nca Nco)) {
        die "Required column $required not found in $source_file header: $header\n"
          unless exists $idx{$required};
    }

    my ($frq_a_col) = grep { /^FRQ_A_/ } @cols;
    my ($frq_u_col) = grep { /^FRQ_U_/ } @cols;
    my $neff_col = exists $idx{Neff} ? 'Neff' : '';
    my $neff_half_col = exists $idx{Neff_half} ? 'Neff_half' : '';

    my $rows = 0;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^\s*$/;

        my @v = split /\t/, $line, -1;
        @v = split /\s+/, $line if @v == 1;

        my $chr = value(\@v, \%idx, 'CHR');
        my $is_chrx = ($chr =~ /^(?:X|23|chrX)$/i) ? 1 : 0;

        my %row = (
            CHR          => normalize_chr($chr),
            BP           => value(\@v, \%idx, 'BP'),
            A1           => value(\@v, \%idx, 'A1'),
            A2           => value(\@v, \%idx, 'A2'),
            SNP          => value(\@v, \%idx, 'SNP'),
            GWAS_TAG     => $tag,
            SOURCE_FILE  => $source_file,
            CHR_ORIGINAL => $chr,
            IS_CHRX      => $is_chrx,
            FRQ_A        => $frq_a_col ? value(\@v, \%idx, $frq_a_col) : '',
            FRQ_U        => $frq_u_col ? value(\@v, \%idx, $frq_u_col) : '',
            INFO         => value(\@v, \%idx, 'INFO'),
            OR           => value(\@v, \%idx, 'OR'),
            SE           => value(\@v, \%idx, 'SE'),
            P            => value(\@v, \%idx, 'P'),
            ngt          => value(\@v, \%idx, 'ngt'),
            Direction    => value(\@v, \%idx, 'Direction'),
            HetISqt      => value(\@v, \%idx, 'HetISqt'),
            HetDf        => value(\@v, \%idx, 'HetDf'),
            HetPVa       => value(\@v, \%idx, 'HetPVa'),
            Nca          => value(\@v, \%idx, 'Nca'),
            Nco          => value(\@v, \%idx, 'Nco'),
            Neff         => $neff_col ? value(\@v, \%idx, $neff_col) : '',
            Neff_half    => $neff_half_col ? value(\@v, \%idx, $neff_half_col) : '',
        );

        print {$out} join("\t", map { defined $row{$_} ? $row{$_} : '' } @out_cols), "\n";
        $rows++;
        last if $limit && $rows >= $limit;
    }

    close $fh;
    return { rows => $rows, header => $header, status => 'OK' };
}

sub value {
    my ($vals, $idx, $col) = @_;
    return '' unless exists $idx->{$col};
    return $vals->[ $idx->{$col} ] // '';
}

sub normalize_chr {
    my ($chr) = @_;
    $chr =~ s/^chr//i;
    return 'X' if $chr =~ /^(?:23|X)$/i;
    return $chr;
}
