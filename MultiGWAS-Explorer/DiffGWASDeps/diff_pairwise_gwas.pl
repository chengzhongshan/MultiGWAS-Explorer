#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);
use IO::Compress::Gzip qw($GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);
use JSON::PP qw(decode_json);

my $config_file = '';
my $input = '';
my $output = '';
my $manifest = '';
my $group_tag_col = 'GWAS_TAG';
my $beta_col = 'BETA';
my $se_col = 'SE';
my $p_col = 'P';
my $rho = 0;
my $limit = 0;

GetOptions(
    'config=s'        => \$config_file,
    'input=s'         => \$input,
    'output=s'        => \$output,
    'manifest=s'      => \$manifest,
    'group-tag-col=s' => \$group_tag_col,
    'beta-col=s'      => \$beta_col,
    'se-col=s'        => \$se_col,
    'p-col=s'         => \$p_col,
    'rho=f'           => \$rho,
    'limit=i'         => \$limit,
) or die usage();

die "--config is required\n" unless length $config_file;
my $cfg = load_json($config_file);

$input = cfg_or($cfg, 'input', $input);
$output = cfg_or($cfg, 'output', $output);
$manifest = cfg_or($cfg, 'manifest', $manifest);
$group_tag_col = cfg_or($cfg, 'group_tag_col', $group_tag_col);
$beta_col = cfg_or($cfg, 'beta_col', $beta_col);
$se_col = cfg_or($cfg, 'se_col', $se_col);
$p_col = cfg_or($cfg, 'p_col', $p_col);
$rho = cfg_or($cfg, 'rho', $rho);
$limit = cfg_or($cfg, 'limit', $limit);

die "Config $config_file must define input\n" unless length $input;
die "Config $config_file must define output\n" unless length $output;
die "Config $config_file must define manifest\n" unless length $manifest;
die "rho must be greater than -1 and less than 1\n" unless $rho > -1 && $rho < 1;
die "Config $config_file must define pairs as a non-empty object\n"
  unless ref($cfg->{pairs}) eq 'HASH' && keys %{ $cfg->{pairs} };

my @base_cols = @{ $cfg->{base_cols} && ref($cfg->{base_cols}) eq 'ARRAY' ? $cfg->{base_cols} : [qw(CHR BP A1 A2 SNP)] };

my @out_cols = qw(
  CHR BP A1 A2 SNP PAIR_TAG
  GROUP1_GWAS_TAG GROUP2_GWAS_TAG
  GROUP1_SOURCE_FILE GROUP2_SOURCE_FILE
  GROUP1_BETA GROUP2_BETA DIFF_BETA
  GROUP1_SE GROUP2_SE DIFF_SE
  GROUP1_Z GROUP2_Z DIFF_Z DIFF_P
  GROUP1_P GROUP2_P
  GROUP1_FRQ_A GROUP1_FRQ_U GROUP2_FRQ_A GROUP2_FRQ_U
  GROUP1_INFO GROUP2_INFO
  CHR_ORIGINAL IS_CHRX
);

my %pair_defs;
for my $pair_tag (sort keys %{ $cfg->{pairs} }) {
    my $def = $cfg->{pairs}{$pair_tag};
    die "Pair $pair_tag must be a 2-element array\n"
      unless ref($def) eq 'ARRAY' && @{$def} == 2;
    $pair_defs{$pair_tag} = [ @{$def} ];
}

my ($header, $idx_ref) = read_header($input);
my %idx = %{$idx_ref};
for my $required (@base_cols, $group_tag_col, $beta_col, $se_col, $p_col, qw(SOURCE_FILE CHR_ORIGINAL IS_CHRX)) {
    die "Required column $required not found in header\n" unless exists $idx{$required};
}

my %stats = (
    rows_read          => 0,
    groups_seen        => 0,
    pairs_written      => 0,
    skipped_no_pair    => 0,
    skipped_bad_num    => 0,
    skipped_bad_var    => 0,
    skipped_unknown_gw => 0,
);

my $in = open_reader($input);
my $out = open_writer($output);
print {$out} join("\t", @out_cols), "\n";
<$in>;

my $current_key = '';
my @bucket;
while (my $line = <$in>) {
    chomp $line;
    $line =~ s/\r$//;
    next if $line =~ /^\s*$/;
    my @v = split /\t/, $line, -1;
    my $key = join("\t", map { $v[ $idx{$_} ] // '' } @base_cols);
    if (@bucket && $key ne $current_key) {
        process_bucket(\@bucket);
        @bucket = ();
        last if $limit && $stats{pairs_written} >= $limit;
    }
    $current_key = $key;
    push @bucket, \@v;
    $stats{rows_read}++;
}
process_bucket(\@bucket) if @bucket && (!$limit || $stats{pairs_written} < $limit);

close $in or die "Failed closing input $input: $!\n";
close $out or die "Failed closing output $output: $!\n";

open my $man, '>', $manifest or die "Cannot write $manifest: $!\n";
print {$man} join("\t", qw(METRIC VALUE)), "\n";
for my $metric (sort keys %stats) {
    print {$man} join("\t", $metric, $stats{$metric}), "\n";
}
print {$man} join("\t", 'input', $input), "\n";
print {$man} join("\t", 'output', $output), "\n";
print {$man} join("\t", 'rho', $rho), "\n";
print {$man} join("\t", 'pairs', join(',', map { $_ . '=' . join(':', @{ $pair_defs{$_} }) } sort keys %pair_defs)), "\n";
close $man;

warn "Differential GWAS output: $output\n";
warn "Manifest: $manifest\n";
warn "Pairs written: $stats{pairs_written}\n";

sub process_bucket {
    my ($bucket) = @_;
    return unless @{$bucket};
    $stats{groups_seen}++;

    my %rows_by_group;
    for my $row (@{$bucket}) {
        my $tag = value($row, $group_tag_col);
        if (!defined $tag || $tag eq '') {
            $stats{skipped_unknown_gw}++;
            next;
        }
        $rows_by_group{$tag} ||= $row;
    }

    for my $pair_tag (sort keys %pair_defs) {
        my ($g1_tag, $g2_tag) = @{ $pair_defs{$pair_tag} };
        my $g1 = $rows_by_group{$g1_tag};
        my $g2 = $rows_by_group{$g2_tag};

        if (!$g1 || !$g2) {
            $stats{skipped_no_pair}++;
            next;
        }

        my $b1 = numeric(value($g1, $beta_col));
        my $b2 = numeric(value($g2, $beta_col));
        my $s1 = numeric(value($g1, $se_col));
        my $s2 = numeric(value($g2, $se_col));
        if (!defined $b1 || !defined $b2 || !defined $s1 || !defined $s2 || $s1 <= 0 || $s2 <= 0) {
            $stats{skipped_bad_num}++;
            next;
        }

        my $diff_beta = $b1 - $b2;
        my $var = $s1 * $s1 + $s2 * $s2 - 2 * $rho * $s1 * $s2;
        if ($var <= 0) {
            $stats{skipped_bad_var}++;
            next;
        }

        my $diff_se = sqrt($var);
        my $z1 = $b1 / $s1;
        my $z2 = $b2 / $s2;
        my $diff_z = $diff_beta / $diff_se;
        my $diff_p = two_sided_p_from_z($diff_z);

        my %o = (
            CHR               => sas_chr(value($g1, 'CHR')),
            BP                => value($g1, 'BP'),
            A1                => value($g1, 'A1'),
            A2                => value($g1, 'A2'),
            SNP               => value($g1, 'SNP'),
            PAIR_TAG          => $pair_tag,
            GROUP1_GWAS_TAG   => $g1_tag,
            GROUP2_GWAS_TAG   => $g2_tag,
            GROUP1_SOURCE_FILE => value($g1, 'SOURCE_FILE'),
            GROUP2_SOURCE_FILE => value($g2, 'SOURCE_FILE'),
            GROUP1_BETA       => fmt($b1),
            GROUP2_BETA       => fmt($b2),
            DIFF_BETA         => fmt($diff_beta),
            GROUP1_SE         => fmt($s1),
            GROUP2_SE         => fmt($s2),
            DIFF_SE           => fmt($diff_se),
            GROUP1_Z          => fmt($z1),
            GROUP2_Z          => fmt($z2),
            DIFF_Z            => fmt($diff_z),
            DIFF_P            => p_fmt($diff_p),
            GROUP1_P          => value($g1, $p_col),
            GROUP2_P          => value($g2, $p_col),
            GROUP1_FRQ_A      => value($g1, 'FRQ_A'),
            GROUP1_FRQ_U      => value($g1, 'FRQ_U'),
            GROUP2_FRQ_A      => value($g2, 'FRQ_A'),
            GROUP2_FRQ_U      => value($g2, 'FRQ_U'),
            GROUP1_INFO       => value($g1, 'INFO'),
            GROUP2_INFO       => value($g2, 'INFO'),
            CHR_ORIGINAL      => value($g1, 'CHR_ORIGINAL'),
            IS_CHRX           => value($g1, 'IS_CHRX'),
        );

        print {$out} join("\t", map { defined $o{$_} ? $o{$_} : '' } @out_cols), "\n";
        $stats{pairs_written}++;
        last if $limit && $stats{pairs_written} >= $limit;
    }
}

sub read_header {
    my ($path) = @_;
    my $fh = open_reader($path);
    my $h = <$fh>;
    close $fh;
    die "Input is empty: $path\n" unless defined $h;
    chomp $h;
    $h =~ s/\r$//;
    $h =~ s/^#//;
    my @cols = split /\t/, $h, -1;
    my %idx;
    for my $i (0 .. $#cols) {
        $idx{$cols[$i]} = $i;
    }
    return ($h, \%idx);
}

sub open_reader {
    my ($path) = @_;
    my $fh;
    if ($path =~ /\.gz$/i) {
        if (prefer_shell_gzip($path)) {
            my $cmd = 'gzip -dc ' . shell_quote($path);
            open $fh, '-|', $cmd or die "Cannot read gzip input $path with gzip -dc: $!\n";
        }
        else {
            $fh = IO::Uncompress::Gunzip->new($path)
              or die "Cannot open gzip input $path: $GunzipError\n";
        }
    }
    else {
        open $fh, '<', $path or die "Cannot read $path: $!\n";
    }
    return $fh;
}

sub open_writer {
    my ($path) = @_;
    my $fh;
    if ($path =~ /\.gz$/i) {
        if (prefer_shell_gzip($path)) {
            my $cmd = 'gzip -c > ' . shell_quote($path);
            open $fh, '|-', $cmd or die "Cannot write gzip output $path with gzip -c: $!\n";
        }
        else {
            $fh = IO::Compress::Gzip->new($path)
              or die "Cannot write gzip output $path: $GzipError\n";
        }
    }
    else {
        open $fh, '>', $path or die "Cannot write $path: $!\n";
    }
    return $fh;
}

sub value {
    my ($row, $col) = @_;
    return '' unless exists $idx{$col};
    return $row->[ $idx{$col} ] // '';
}

sub numeric {
    my ($x) = @_;
    return undef unless defined $x;
    return undef if $x eq '' || $x =~ /^(?:NA|NaN|null|\.)$/i;
    return undef unless $x =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
    return 0 + $x;
}

sub two_sided_p_from_z {
    my ($z) = @_;
    return erfc_approx(abs($z) / sqrt(2));
}

sub erfc_approx {
    my ($x) = @_;
    my $z = abs($x);
    my $t = 1 / (1 + 0.5 * $z);
    my $poly = ((((((((0.17087277 * $t - 0.82215223) * $t + 1.48851587) * $t
      - 1.13520398) * $t + 0.27886807) * $t - 0.18628806) * $t
      + 0.09678418) * $t + 0.37409196) * $t + 1.00002368) * $t;
    my $r = $t * exp(-$z * $z - 1.26551223 + $poly);
    return $x >= 0 ? $r : 2 - $r;
}

sub fmt {
    my ($x) = @_;
    return '' unless defined $x;
    return sprintf('%.10g', $x);
}

sub p_fmt {
    my ($x) = @_;
    return '' unless defined $x;
    return sprintf('%.6e', $x);
}

sub sas_chr {
    my ($chr) = @_;
    return '' unless defined $chr;
    return '23' if $chr =~ /^(?:X|chrX)$/i;
    return $chr;
}

sub cfg_or {
    my ($cfg, $key, $fallback) = @_;
    return $fallback unless exists $cfg->{$key};
    return $cfg->{$key};
}

sub prefer_shell_gzip {
    my ($path) = @_;
    return 1 if defined $path && $path =~ m{^/mnt/};
    return 0;
}

sub shell_quote {
    my ($text) = @_;
    return "''" unless defined $text && length $text;
    $text =~ s/'/'\"'\"'/g;
    return "'$text'";
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
  perl diff_pairwise_gwas.pl --config diff_config.json [options]

Options:
  --config FILE.json    Required JSON config with input, output, manifest, pairs
  --input FILE.tsv.gz   Override config input
  --output FILE.tsv.gz  Override config output
  --manifest FILE.tsv   Override config manifest
  --group-tag-col NAME  Column containing merged long GWAS tag. Default: GWAS_TAG
  --beta-col NAME       Effect-size column in merged long input. Default: BETA
  --se-col NAME         SE column in merged long input. Default: SE
  --p-col NAME          P-value column in merged long input. Default: P
  --rho FLOAT           Correlation between paired betas. Default: 0
  --limit N             Optional output-pair limit for debugging
USAGE
}
