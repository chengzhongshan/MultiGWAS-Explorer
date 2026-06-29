#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
BEGIN {
    my $self_dir = __FILE__;
    $self_dir =~ s{\\}{/}g;
    $self_dir =~ s{/[^/]+$}{};
    $self_dir = '.' unless length $self_dir;
    require lib;
    lib->import($self_dir);
}
use Getopt::Long qw(GetOptions);
use IO::Compress::Gzip qw($GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);
use JSON::PP qw(decode_json);
use POSIX qw(log);
use DiffGWASRawSchema qw(
  normalize_header_name
  resolve_raw_header_aliases
);

my $EXTERNAL_GZIP = '';

my $config_file = '';
my $input_dir = '';
my $output = '';
my $manifest = '';
my $limit = 0;

GetOptions(
    'config=s'    => \$config_file,
    'input-dir=s' => \$input_dir,
    'output=s'    => \$output,
    'manifest=s'  => \$manifest,
    'limit=i'     => \$limit,
) or die usage();

die "--config is required\n" unless length $config_file;
my $cfg = load_json($config_file);

$input_dir = cfg_or($cfg, 'input_dir', $input_dir);
$output = cfg_or($cfg, 'output', $output);
$manifest = cfg_or($cfg, 'manifest', $manifest);
$limit = cfg_or($cfg, 'limit', $limit);

die "Config $config_file must define input_dir\n" unless length $input_dir;
die "Config $config_file must define output\n" unless length $output;
die "Config $config_file must define manifest\n" unless length $manifest;
die "Config $config_file must define groups as a non-empty array\n"
  unless ref($cfg->{groups}) eq 'ARRAY' && @{ $cfg->{groups} };
my $raw_column_aliases = cfg_or($cfg, 'raw_column_aliases', {});

my @out_cols = qw(
  CHR BP A1 A2 SNP GWAS_TAG SOURCE_FILE CHR_ORIGINAL IS_CHRX
  FRQ_A FRQ_U INFO BETA SE P NCAS NCON NEFF
);

my $out = open_writer($output);
print {$out} join("\t", @out_cols), "\n";

open my $man, '>', $manifest or die "Cannot write $manifest: $!\n";
print {$man} join("\t", qw(GWAS_TAG SOURCE_FILE ROWS_WRITTEN HEADER STATUS)), "\n";

for my $group (@{ $cfg->{groups} }) {
    die "Each group must be a JSON object\n" unless ref($group) eq 'HASH';
    my $tag = $group->{tag} // '';
    die "Group tag is required in $config_file\n" unless length $tag;
    my $files = $group->{files};
    die "Group $tag must define a non-empty files array\n"
      unless ref($files) eq 'ARRAY' && @{$files};

    for my $file (@{$files}) {
        my $path = resolve_path($input_dir, $file);
        unless (-e $path) {
            print {$man} join("\t", $tag, $file, 0, '', 'MISSING'), "\n";
            warn "Missing expected file: $path\n";
            next;
        }

        my $rows = process_file($path, $file, $tag, $out, $limit);
        print {$man} join("\t", $tag, $file, $rows->{rows}, $rows->{header}, $rows->{status}), "\n";
        warn "Finished $file as $tag: $rows->{rows} rows\n";
    }
}

close $out or die "Failed closing output $output: $!\n";
close $man or die "Failed closing manifest $manifest: $!\n";

warn "Merged long-format table: $output\n";
warn "Manifest: $manifest\n";

sub process_file {
    my ($path, $source_file, $tag, $out_fh, $limit_rows) = @_;

    my $fh = open_reader($path);

    my ($header, @cols);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^##/;
        $header = $line;
        @cols = split /\t/, $header, -1;
        last;
    }

    unless (defined $header) {
        close $fh;
        return { rows => 0, header => '', status => 'EMPTY_OR_UNREADABLE' };
    }

    my %idx;
    for my $i (0 .. $#cols) {
        $idx{$cols[$i]} = $i;
    }
    my %resolved = resolve_header_aliases(\@cols, \%idx, $source_file, $header, $raw_column_aliases);

    my $rows = 0;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^\s*$/;

        my @v = split /\t/, $line, -1;
        my $chr = value(\@v, \%resolved, 'CHROM');
        my $norm_chr = normalize_chr($chr);
        my $is_chrx = ($norm_chr =~ /^(?:X|23)$/i) ? 1 : 0;

        my %row = (
            CHR          => $norm_chr,
            BP           => value(\@v, \%resolved, 'POS'),
            A1           => value(\@v, \%resolved, 'A1'),
            A2           => value(\@v, \%resolved, 'A2'),
            SNP          => value(\@v, \%resolved, 'ID'),
            GWAS_TAG     => $tag,
            SOURCE_FILE  => $source_file,
            CHR_ORIGINAL => $chr,
            IS_CHRX      => $is_chrx,
            FRQ_A        => value(\@v, \%resolved, 'FCAS'),
            FRQ_U        => value(\@v, \%resolved, 'FCON'),
            INFO         => value(\@v, \%resolved, 'IMPINFO'),
            BETA         => effect_value(\@v, \%resolved, \@cols),
            SE           => value(\@v, \%resolved, 'SE'),
            P            => value(\@v, \%resolved, 'PVAL'),
            NCAS         => value(\@v, \%resolved, 'NCAS'),
            NCON         => value(\@v, \%resolved, 'NCON'),
            NEFF         => value(\@v, \%resolved, 'NEFF'),
        );

        print {$out_fh} join("\t", map { defined $row{$_} ? $row{$_} : '' } @out_cols), "\n";
        $rows++;
        last if $limit_rows && $rows >= $limit_rows;
    }

    close $fh;
    return { rows => $rows, header => $header, status => 'OK' };
}

sub open_reader {
    my ($path) = @_;
    my $fh;
    if ($path =~ /\.gz$/i) {
        my $src = $path =~ /^[A-Za-z]:[\\\/]/ ? $path : cygpath_to_win($path);
        my $gzip = external_gzip_path();
        if ($gzip && open($fh, '-|', $gzip, '-dc', $src)) {
            return $fh;
        }
        $fh = IO::Uncompress::Gunzip->new($src)
          or die "Cannot open gzip input $src: $GunzipError\n";
    }
    else {
        my $src = $path =~ /^[A-Za-z]:[\\\/]/ ? $path : cygpath_to_win($path);
        open $fh, '<', $src or die "Cannot read $src: $!\n";
    }
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

sub open_writer {
    my ($path) = @_;
    my $fh;
    if ($path =~ /\.gz$/i) {
        $fh = IO::Compress::Gzip->new($path, Level => 1)
          or die "Cannot write gzip output $path: $GzipError\n";
    }
    else {
        open $fh, '>', $path or die "Cannot write $path: $!\n";
    }
    return $fh;
}

sub resolve_path {
    my ($dir, $file) = @_;
    return $file if $file =~ m{^(?:[A-Za-z]:[\\/]|/mnt/|/)};
    $dir =~ s{[\\/]\z}{};
    return "$dir/$file";
}

sub cygpath_to_win {
    my ($path) = @_;
    return '' unless defined $path;
    return $path if $path =~ /^[A-Za-z]:[\\\/]/;
    return $path if $^O !~ /^(?:cygwin|MSWin32)$/i;
    if ($path =~ m{^/mnt/([A-Za-z])/(.*)$}) {
        my ($drive, $rest) = ($1, $2);
        $rest =~ s{/}{\\}g;
        return uc($drive) . ":\\" . $rest;
    }
    my $win = $path;
    $win =~ s{/}{\\}g;
    return $win;
}

sub value {
    my ($vals, $idx, $col) = @_;
    return '' unless exists $idx->{$col};
    return $vals->[ $idx->{$col} ] // '';
}

sub effect_value {
    my ($vals, $resolved, $cols) = @_;
    my $val = value($vals, $resolved, 'BETA');
    return '' unless defined $val && length $val;
    return $val unless exists $resolved->{BETA};
    my $actual = $cols->[ $resolved->{BETA} ] // '';
    my $norm = normalize_header($actual);
    if ($norm eq 'OR' || $norm eq 'ODDSRATIO') {
        return '' unless $val =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?$/ && $val > 0;
        return log($val);
    }
    return $val;
}

sub resolve_header_aliases {
    my ($cols, $idx, $source_file, $header, $alias_overrides) = @_;
    my %resolved = resolve_raw_header_aliases(
        cols            => $cols,
        idx             => $idx,
        source_file     => $source_file,
        header          => $header,
        alias_overrides => $alias_overrides,
    );
    my @required = qw(CHROM ID POS A1 A2 BETA SE PVAL);
    my @optional = qw(FCAS FCON IMPINFO NCAS NCON NEFF);

    my @resolved_notes = map {
        my $actual = '';
        if (exists $resolved{$_}) {
            $actual = $cols->[ $resolved{$_} ];
            "$_=$actual";
        }
        else {
            "$_=NA";
        }
    } (@required, @optional);
    warn "Resolved columns for $source_file: " . join(', ', @resolved_notes) . "\n";

    return %resolved;
}

sub normalize_header {
    return normalize_header_name(@_);
}

sub normalize_chr {
    my ($chr) = @_;
    $chr =~ s/^chr//i if defined $chr;
    return '23' if defined $chr && $chr =~ /^X$/i;
    return $chr;
}

sub cfg_or {
    my ($cfg, $key, $fallback) = @_;
    return $fallback unless exists $cfg->{$key};
    return $cfg->{$key};
}

sub load_json {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read config $path: $!\n";
    local $/;
    my $json = <$fh>;
    close $fh;
    my $cfg = decode_json($json);
    die "Config root must be a JSON object: $path\n" unless ref($cfg) eq 'HASH';
    return $cfg;
}

sub usage {
    return <<"USAGE";
Usage:
  perl merge_pgc_vcf_sumstats_long.pl --config merge_config.json [options]

Options:
  --config FILE.json    Required JSON config with input_dir, output, manifest, groups
  --input-dir DIR       Override config input_dir
  --output FILE.tsv.gz  Override config output
  --manifest FILE.tsv   Override config manifest
  --limit N             Optional per-file row limit for debugging

This helper auto-detects common raw GWAS header synonyms. Examples:
  P column      : PVAL, P, PVALUE, P_VALUE, p-value
  SNP column    : ID, SNP, RSID, MARKERNAME
  position      : POS, BP, POSITION
  chromosome    : CHROM, CHR, CHROMOSOME
  effect allele : A1, EA, EFFECT_ALLELE, ALLELE1
  other allele  : A2, NEA, OTHER_ALLELE, ALLELE2, REF

Core required fields after alias matching:
  chromosome, SNP id, position, A1, A2, beta/effect, SE, P

Optional fields after alias matching:
  case/control frequency, info/imputation score, case/control/effective sample size
USAGE
}
