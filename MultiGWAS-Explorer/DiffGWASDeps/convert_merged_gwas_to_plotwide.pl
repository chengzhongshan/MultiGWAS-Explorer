#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);
use IO::Compress::Gzip qw($GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);
use POSIX qw(erfc);
use lib $Bin;
use DiffGWASConfig qw(load_config_file);

my $config_file = '';

GetOptions(
    'config=s' => \$config_file,
) or die usage();

die "--config is required\n" unless length $config_file;
my $cfg = load_config_file($config_file);

my $input = cfg_or($cfg, 'input', '');
my $output = cfg_or($cfg, 'output', '');
my $manifest = cfg_or($cfg, 'manifest', '');
my $base_map = cfg_or($cfg, 'merged_base_cols', {});
my $group_tracks = cfg_or($cfg, 'merged_group_tracks', []);
my $pairs = cfg_or($cfg, 'pairs', []);
my $extra_tracks = cfg_or($cfg, 'merged_extra_tracks', []);
my $wide_columns = cfg_or($cfg, 'wide_columns', []);

$input = local_path($input);
$output = local_path($output);
$manifest = local_path($manifest);

die "input is required in merged-wide config\n" unless defined $input && length $input;
die "output is required in merged-wide config\n" unless defined $output && length $output;
die "manifest is required in merged-wide config\n" unless defined $manifest && length $manifest;
die "merged_base_cols must be a JSON object\n" unless ref($base_map) eq 'HASH';
die "merged_group_tracks must be a non-empty array\n" unless ref($group_tracks) eq 'ARRAY' && @{$group_tracks} >= 2;
die "pairs must be a non-empty array\n" unless ref($pairs) eq 'ARRAY' && @{$pairs};
die "wide_columns must be a non-empty array\n" unless ref($wide_columns) eq 'ARRAY' && @{$wide_columns};

my ($header, $idx_ref) = read_header($input);
my %idx = %{$idx_ref};

my $chr_col = pick_first_defined($base_map->{chr}, find_header_actual([ keys %idx ], [qw(CHR chromosome chrom)]));
my $bp_col  = pick_first_defined($base_map->{bp},  find_header_actual([ keys %idx ], [qw(BP POS position)]));
my $snp_col = pick_first_defined($base_map->{snp}, find_header_actual([ keys %idx ], [qw(SNP rsid marker id)]));
my $a1_col  = pick_first_defined($base_map->{a1},  find_header_actual([ keys %idx ], [qw(A1 EA EFFECT_ALLELE ALT)]));
my $a2_col  = pick_first_defined($base_map->{a2},  find_header_actual([ keys %idx ], [qw(A2 NEA OTHER_ALLELE REF)]));

for my $required ($chr_col, $bp_col, $snp_col) {
    die "Merged input is missing one of the required base columns CHR/BP/SNP\n"
      unless defined $required && exists $idx{$required};
}

my %group_by_id = map {
    my $id = $_->{id} || '';
    $id => $_
} grep { ref($_) eq 'HASH' && defined($_->{id}) && length($_->{id}) } @{$group_tracks};

for my $pair (@{$pairs}) {
    die "Each pair must be an object\n" unless ref($pair) eq 'HASH';
    my $g1 = $pair->{source_group1} || '';
    my $g2 = $pair->{source_group2} || '';
    my $prefix = $pair->{prefix} || '';
    die "Pair is missing source_group1/source_group2/prefix\n" unless length($g1) && length($g2) && length($prefix);
    die "Pair source_group1 $g1 is not defined in merged_group_tracks\n" unless exists $group_by_id{$g1};
    die "Pair source_group2 $g2 is not defined in merged_group_tracks\n" unless exists $group_by_id{$g2};
    for my $need (qw(beta_col se_col p_col)) {
        my $c1 = $group_by_id{$g1}{$need};
        my $c2 = $group_by_id{$g2}{$need};
        die "Merged group track $g1 is missing $need\n" unless defined $c1 && exists $idx{$c1};
        die "Merged group track $g2 is missing $need\n" unless defined $c2 && exists $idx{$c2};
    }
}

for my $track (@{$extra_tracks}) {
    next unless ref($track) eq 'HASH';
    next unless defined($track->{p_col}) && length($track->{p_col});
    next unless defined($track->{z_col}) && length($track->{z_col});
    die "Merged extra-track P column $track->{p_col} not found in input\n" unless exists $idx{ $track->{p_col} };
    die "Merged extra-track Z column $track->{z_col} not found in input\n" unless exists $idx{ $track->{z_col} };
    if (defined($track->{beta_col}) && length($track->{beta_col})) {
        die "Merged extra-track beta column $track->{beta_col} not found in input\n" unless exists $idx{ $track->{beta_col} };
    }
    if (defined($track->{se_col}) && length($track->{se_col})) {
        die "Merged extra-track se column $track->{se_col} not found in input\n" unless exists $idx{ $track->{se_col} };
    }
}

my $in = open_reader($input);
my $out = open_writer($output);
print {$out} join("\t", @{$wide_columns}), "\n";
<$in>;

my %stats = (
    rows_read    => 0,
    rows_written => 0,
    input        => $input,
    output       => $output,
    manifest     => $manifest,
    pair_count   => scalar(@{$pairs}),
    extra_track_count => scalar(@{$extra_tracks}),
);

while (my $line = <$in>) {
    chomp $line;
    $line =~ s/\r$//;
    next if $line =~ /^\s*$/;
    $stats{rows_read}++;
    my @f = split /\t/, $line, -1;
    my %row = map { $_ => '' } @{$wide_columns};
    $row{CHR} = normalize_chr($f[ $idx{$chr_col} ]);
    $row{BP}  = normalize_numeric_text($f[ $idx{$bp_col} ]);
    $row{SNP} = $f[ $idx{$snp_col} ] // '';
    $row{A1}  = defined $a1_col && exists $idx{$a1_col} ? ($f[ $idx{$a1_col} ] // '') : '';
    $row{A2}  = defined $a2_col && exists $idx{$a2_col} ? ($f[ $idx{$a2_col} ] // '') : '';

    for my $pair (@{$pairs}) {
        my $prefix = $pair->{prefix};
        my $g1 = $group_by_id{ $pair->{source_group1} };
        my $g2 = $group_by_id{ $pair->{source_group2} };
        my $b1 = numeric($f[ $idx{ $g1->{beta_col} } ]);
        my $s1 = numeric($f[ $idx{ $g1->{se_col} } ]);
        my $p1 = numeric($f[ $idx{ $g1->{p_col} } ]);
        my $b2 = numeric($f[ $idx{ $g2->{beta_col} } ]);
        my $s2 = numeric($f[ $idx{ $g2->{se_col} } ]);
        my $p2 = numeric($f[ $idx{ $g2->{p_col} } ]);
        $row{"${prefix}_GROUP1_BETA"} = format_numeric($b1);
        $row{"${prefix}_GROUP1_SE"}   = format_numeric($s1);
        $row{"${prefix}_GROUP1_P"}    = format_numeric($p1);
        $row{"${prefix}_GROUP1_Z"}    = format_numeric(z_from_beta_se($b1, $s1));
        $row{"${prefix}_GROUP2_BETA"} = format_numeric($b2);
        $row{"${prefix}_GROUP2_SE"}   = format_numeric($s2);
        $row{"${prefix}_GROUP2_P"}    = format_numeric($p2);
        $row{"${prefix}_GROUP2_Z"}    = format_numeric(z_from_beta_se($b2, $s2));
        my $diff_beta = (defined($b1) && defined($b2)) ? ($b1 - $b2) : undef;
        my $diff_se = diff_se_independent($s1, $s2);
        my $diff_z = z_from_beta_se($diff_beta, $diff_se);
        $row{"${prefix}_DIFF_BETA"}   = format_numeric($diff_beta);
        $row{"${prefix}_DIFF_SE"}     = format_numeric($diff_se);
        $row{"${prefix}_DIFF_P"}      = format_numeric(two_sided_p_from_z($diff_z));
        $row{"${prefix}_STD_DIFF_Z"}  = format_numeric($diff_z);
        $row{"${prefix}_STD_DIFF_P"}  = format_numeric(two_sided_p_from_z($diff_z));
    }

    for my $track (@{$extra_tracks}) {
        next unless ref($track) eq 'HASH';
        my $id = $track->{id} || next;
        my $p = numeric($f[ $idx{ $track->{p_col} } ]);
        my $z = numeric($f[ $idx{ $track->{z_col} } ]);
        $row{"${id}_P"} = format_numeric($p) if exists $row{"${id}_P"};
        $row{"${id}_Z"} = format_numeric($z) if exists $row{"${id}_Z"};
        if (defined($track->{beta_col}) && length($track->{beta_col}) && exists $row{"${id}_BETA"}) {
            $row{"${id}_BETA"} = format_numeric(numeric($f[ $idx{ $track->{beta_col} } ]));
        }
        if (defined($track->{se_col}) && length($track->{se_col}) && exists $row{"${id}_SE"}) {
            $row{"${id}_SE"} = format_numeric(numeric($f[ $idx{ $track->{se_col} } ]));
        }
    }

    print {$out} join("\t", map { defined $row{$_} ? $row{$_} : '' } @{$wide_columns}), "\n";
    $stats{rows_written}++;
}

close $in or die "Failed closing input $input: $!\n";
close $out or die "Failed closing output $output: $!\n";

open my $mf, '>', $manifest or die "Cannot write $manifest: $!\n";
print {$mf} join("\t", qw(METRIC VALUE)), "\n";
for my $metric (qw(input output manifest rows_read rows_written pair_count extra_track_count)) {
    print {$mf} join("\t", $metric, $stats{$metric}), "\n";
}
print {$mf} join("\t", 'columns', join(',', @{$wide_columns})), "\n";
close $mf or die "Cannot close $manifest: $!\n";

print "Input:        $input\n";
print "Output:       $output\n";
print "Manifest:     $manifest\n";
print "Rows read:    $stats{rows_read}\n";
print "Rows written: $stats{rows_written}\n";

sub cfg_or {
    my ($cfg, $key, $fallback) = @_;
    return $fallback unless ref($cfg) eq 'HASH' && exists $cfg->{$key};
    return $cfg->{$key};
}

sub open_reader {
    my ($path) = @_;
    $path = local_path($path);
    if ($path =~ /\.(?:gz|bgz|bgzip)$/i) {
        my $fh = IO::Uncompress::Gunzip->new($path)
          or die "Cannot open gzip input $path: $GunzipError\n";
        return $fh;
    }
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    return $fh;
}

sub open_writer {
    my ($path) = @_;
    $path = local_path($path);
    my $dir = dirname($path);
    make_path($dir) if defined $dir && length $dir && !-d $dir;
    if ($path =~ /\.(?:gz|bgz|bgzip)$/i) {
        my $fh = IO::Compress::Gzip->new($path)
          or die "Cannot write gzip output $path: $GzipError\n";
        return $fh;
    }
    open my $fh, '>', $path or die "Cannot write $path: $!\n";
    return $fh;
}

sub read_header {
    my ($path) = @_;
    my $fh = open_reader($path);
    my $header = <$fh>;
    die "Input file is empty: $path\n" unless defined $header;
    chomp $header;
    $header =~ s/\r$//;
    my @cols = split /\t/, $header, -1;
    close $fh;
    my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
    return ($header, \%idx);
}

sub local_path {
    my ($path) = @_;
    return '' unless defined $path;
    if ($^O =~ /^(?:cygwin|MSWin32)$/i) {
        if ($path =~ m{^/mnt/([A-Za-z])/(.*)$}) {
            my ($drive, $rest) = ($1, $2);
            $rest =~ s{/}{\\}g;
            return uc($drive) . ":\\" . $rest;
        }
        if ($path =~ m{^/cygdrive/([A-Za-z])/(.*)$}) {
            my ($drive, $rest) = ($1, $2);
            $rest =~ s{/}{\\}g;
            return uc($drive) . ":\\" . $rest;
        }
    }
    return $path;
}

sub find_header_actual {
    my ($cols, $aliases) = @_;
    my %lookup = map { uc($_) => $_ } @{$cols || []};
    for my $alias (@{$aliases || []}) {
        my $u = uc($alias // '');
        return $lookup{$u} if exists $lookup{$u};
    }
    return undef;
}

sub pick_first_defined {
    for my $v (@_) {
        return $v if defined $v && $v ne '';
    }
    return undef;
}

sub numeric {
    my ($value) = @_;
    return undef unless defined $value;
    $value =~ s/^\s+|\s+$//g;
    return undef unless length $value;
    return undef if $value =~ /^(?:NA|N\/A|null|\.)$/i;
    return $value + 0 if $value =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?$/;
    return undef;
}

sub normalize_numeric_text {
    my ($value) = @_;
    my $num = numeric($value);
    return '' unless defined $num;
    return format_numeric($num);
}

sub format_numeric {
    my ($value) = @_;
    return '' unless defined $value;
    return sprintf('%.12g', $value);
}

sub normalize_chr {
    my ($value) = @_;
    $value = '' unless defined $value;
    $value =~ s/^\s+|\s+$//g;
    $value =~ s/^chr//i;
    return $value;
}

sub z_from_beta_se {
    my ($beta, $se) = @_;
    return undef unless defined $beta && defined $se && $se > 0;
    return $beta / $se;
}

sub diff_se_independent {
    my ($se1, $se2) = @_;
    return undef unless defined $se1 && defined $se2 && $se1 >= 0 && $se2 >= 0;
    return sqrt(($se1 * $se1) + ($se2 * $se2));
}

sub two_sided_p_from_z {
    my ($z) = @_;
    return undef unless defined $z;
    my $az = abs($z) / sqrt(2);
    my $p = erfc($az);
    return undef unless defined $p;
    return $p;
}

sub usage {
    return <<"USAGE";
Usage:
  perl convert_merged_gwas_to_plotwide.pl --config merged_plotwide.json

The config should describe:
  - input / output / manifest
  - merged_base_cols
  - merged_group_tracks
  - pairs
  - merged_extra_tracks
  - wide_columns
USAGE
}
