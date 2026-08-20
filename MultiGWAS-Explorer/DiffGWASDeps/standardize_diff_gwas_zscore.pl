#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);

my $input =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.tsv.gz';
my $output =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.tsv.gz';
my $manifest =
  '/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.manifest.tsv';
my $z_col = 'DIFF_Z';
my $htsbin = "$Bin/../local/bin";
my $index_output = 1;
my $method = 'mean_sd';
my $clip_lower_quantile = 0.001;
my $clip_upper_quantile = 0.999;

GetOptions(
    'input=s'               => \$input,
    'output=s'              => \$output,
    'manifest=s'            => \$manifest,
    'z-col=s'               => \$z_col,
    'htsbin=s'              => \$htsbin,
    'method=s'              => \$method,
    'clip-lower-quantile=s' => \$clip_lower_quantile,
    'clip-upper-quantile=s' => \$clip_upper_quantile,
    'index-output!'         => \$index_output,
) or die usage();

die "Input file not found: $input\n" unless -s $input;
validate_method($method);
$clip_lower_quantile = validate_quantile($clip_lower_quantile, 'clip-lower-quantile');
$clip_upper_quantile = validate_quantile($clip_upper_quantile, 'clip-upper-quantile');
die "--clip-lower-quantile must be smaller than --clip-upper-quantile\n"
    unless $clip_lower_quantile < $clip_upper_quantile;

my ($header, $cols, $idx) = read_header($input);
my %idx = %$idx;
die "Column $z_col not found in $input\n" unless exists $idx{$z_col};
my $bgzip = resolve_hts_tool($htsbin, 'bgzip');
my $tabix = resolve_hts_tool($htsbin, 'tabix');
my $can_index_output = (
    $index_output
    && defined $bgzip
    && defined $tabix
    && $output =~ /\.gz$/i
    && exists $idx{CHR}
    && exists $idx{BP}
);

my $stats = compute_standardization_stats(
    path                => $input,
    col_i               => $idx{$z_col},
    method              => $method,
    clip_lower_quantile => $clip_lower_quantile,
    clip_upper_quantile => $clip_upper_quantile,
);
my ($n, $mean, $sd) = @{$stats}{qw(n mean sd)};
die "No numeric values found in $z_col\n" unless $n > 0;
die "Cannot standardize because standard deviation is zero\n" unless $sd > 0;

open my $in,  '-|', "zcat '$input'"       or die "Cannot read $input with zcat: $!\n";
my $writer_cmd = $can_index_output
  ? "$bgzip -c > '$output'"
  : "gzip -c > '$output'";
open my $out, '|-', $writer_cmd or die "Cannot write $output: $!\n";

my $raw_header = <$in>;
chomp $raw_header;
$raw_header =~ s/\r$//;
print {$out} join("\t", $raw_header, 'ORIG_DIFF_Z', 'ORIG_DIFF_P', 'STD_DIFF_Z', 'STD_DIFF_P'), "\n";

my $rows = 0;
my $numeric_rows = 0;
while (my $line = <$in>) {
    chomp $line;
    $line =~ s/\r$//;
    my @v = split /\t/, $line, -1;
    my $z = numeric($v[ $idx{$z_col} ]);
    if (defined $z) {
        my $std_z = ($z - $mean) / $sd;
        my $std_p = two_sided_p_from_z($std_z);
        print {$out} join("\t", @v, fmt($z), value(\@v, \%idx, 'DIFF_P'), fmt($std_z), p_fmt($std_p)), "\n";
        $numeric_rows++;
    }
    else {
        print {$out} join("\t", @v, '', value(\@v, \%idx, 'DIFF_P'), '', ''), "\n";
    }
    $rows++;
}

close $in;
close $out or die "Failed closing gzip output $output: $!\n";

my $index_status = 'disabled';
if ($can_index_output) {
    my $seq_col = $idx{CHR} + 1;
    my $bp_col  = $idx{BP} + 1;
    if (system($tabix, '-f', '-s', $seq_col, '-b', $bp_col, '-e', $bp_col, '-S', 1, $output) == 0) {
        $index_status = 'created';
    }
    else {
        warn "Warning: tabix failed for $output; continuing without an index.\n";
        unlink "$output.tbi" if -e "$output.tbi";
        unlink "$output.csi" if -e "$output.csi";
        $index_status = 'failed';
    }
}
elsif ($index_output) {
    $index_status = 'unavailable';
}

open my $man, '>', $manifest or die "Cannot write $manifest: $!\n";
print {$man} join("\t", qw(METRIC VALUE)), "\n";
print {$man} join("\t", 'input',        $input), "\n";
print {$man} join("\t", 'output',       $output), "\n";
print {$man} join("\t", 'z_col',        $z_col), "\n";
print {$man} join("\t", 'inference_p_column', 'DIFF_P'), "\n";
print {$man} join("\t", 'std_diff_p_use', 'legacy_visualization_only_not_inferential'), "\n";
print {$man} join("\t", 'method',       $method), "\n";
print {$man} join("\t", 'clip_lower_quantile', fmt($clip_lower_quantile)), "\n";
print {$man} join("\t", 'clip_upper_quantile', fmt($clip_upper_quantile)), "\n";
print {$man} join("\t", 'clip_lower_value', fmt($stats->{clip_lower_value})), "\n";
print {$man} join("\t", 'clip_upper_value', fmt($stats->{clip_upper_value})), "\n";
print {$man} join("\t", 'lower_tail_clipped_n', $stats->{lower_tail_clipped_n} // 0), "\n";
print {$man} join("\t", 'upper_tail_clipped_n', $stats->{upper_tail_clipped_n} // 0), "\n";
print {$man} join("\t", 'index_output', $index_output ? 1 : 0), "\n";
print {$man} join("\t", 'index_status', $index_status), "\n";
print {$man} join("\t", 'numeric_n',    $n), "\n";
print {$man} join("\t", 'mean',         fmt($mean)), "\n";
print {$man} join("\t", 'sample_sd',    fmt($sd)), "\n";
print {$man} join("\t", 'rows_written', $rows), "\n";
print {$man} join("\t", 'numeric_rows', $numeric_rows), "\n";
close $man;

print "Input:        $input\n";
print "Output:       $output\n";
print "Manifest:     $manifest\n";
print "Method:       $method\n";
print "Index status: $index_status\n";
print "Numeric N:    $n\n";
print "Mean $z_col:  ", fmt($mean), "\n";
print "SD $z_col:    ", fmt($sd), "\n";
if ($method eq 'mean_sd_clipped') {
    print "Clip lower q: ", fmt($clip_lower_quantile), " => ", fmt($stats->{clip_lower_value}), "\n";
    print "Clip upper q: ", fmt($clip_upper_quantile), " => ", fmt($stats->{clip_upper_value}), "\n";
    print "Lower clipped:", ' ', ($stats->{lower_tail_clipped_n} // 0), "\n";
    print "Upper clipped:", ' ', ($stats->{upper_tail_clipped_n} // 0), "\n";
}
print "Rows written: $rows\n";
print "Index:        $output.tbi\n" if $index_status eq 'created';

sub resolve_hts_tool {
    my ($dir, $tool) = @_;
    for my $candidate (
        (defined $dir && length $dir ? ("$dir/$tool", "$dir/$tool.exe") : ()),
        $tool,
    ) {
        next unless defined $candidate && length $candidate;
        return $candidate if -x $candidate || command_exists($candidate);
    }
    return undef;
}

sub command_exists {
    my ($cmd) = @_;
    return 0 unless defined $cmd && length $cmd;
    return scalar(`command -v '$cmd' 2>/dev/null`) ? 1 : 0;
}

sub read_header {
    my ($path) = @_;
    open my $fh, '-|', "zcat '$path'" or die "Cannot read $path with zcat: $!\n";
    my $h = <$fh>;
    close $fh;
    die "Input is empty: $path\n" unless defined $h;
    chomp $h;
    $h =~ s/\r$//;
    $h =~ s/^#//;
    my @c = split /\t/, $h, -1;
    my %i;
    for my $j (0 .. $#c) {
        $i{$c[$j]} = $j;
    }
    return ($h, \@c, \%i);
}

sub compute_standardization_stats {
    my (%args) = @_;
    my $method = $args{method} // 'mean_sd';
    if ($method eq 'mean_sd') {
        my ($n, $mean, $sd) = z_stats($args{path}, $args{col_i});
        return {
            n                   => $n,
            mean                => $mean,
            sd                  => $sd,
            clip_lower_value    => undef,
            clip_upper_value    => undef,
            lower_tail_clipped_n => 0,
            upper_tail_clipped_n => 0,
        };
    }
    if ($method eq 'mean_sd_clipped') {
        return z_stats_clipped(
            $args{path},
            $args{col_i},
            $args{clip_lower_quantile},
            $args{clip_upper_quantile},
        );
    }
    die "Unsupported standardization method: $method\n";
}

sub z_stats {
    my ($path, $col_i) = @_;
    open my $fh, '-|', "zcat '$path'" or die "Cannot read $path with zcat: $!\n";
    <$fh>;
    my ($n, $mean, $m2) = (0, 0, 0);
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        my @v = split /\t/, $line, -1;
        my $x = numeric($v[$col_i]);
        next unless defined $x;
        $n++;
        my $delta = $x - $mean;
        $mean += $delta / $n;
        my $delta2 = $x - $mean;
        $m2 += $delta * $delta2;
    }
    close $fh;
    my $sd = $n > 1 ? sqrt($m2 / ($n - 1)) : 0;
    return ($n, $mean, $sd);
}

sub z_stats_clipped {
    my ($path, $col_i, $lower_q, $upper_q) = @_;
    my @values = read_numeric_values($path, $col_i);
    my $n = scalar @values;
    return {
        n                    => 0,
        mean                 => 0,
        sd                   => 0,
        clip_lower_value     => undef,
        clip_upper_value     => undef,
        lower_tail_clipped_n => 0,
        upper_tail_clipped_n => 0,
    } unless $n;

    @values = sort { $a <=> $b } @values;
    my $lower_value = empirical_quantile(\@values, $lower_q);
    my $upper_value = empirical_quantile(\@values, $upper_q);
    my ($mean, $m2) = (0, 0);
    my ($lower_tail_clipped_n, $upper_tail_clipped_n) = (0, 0);
    my $seen = 0;
    for my $x (@values) {
        my $clipped = $x;
        if ($clipped < $lower_value) {
            $clipped = $lower_value;
            $lower_tail_clipped_n++;
        }
        elsif ($clipped > $upper_value) {
            $clipped = $upper_value;
            $upper_tail_clipped_n++;
        }
        $seen++;
        my $delta = $clipped - $mean;
        $mean += $delta / $seen;
        my $delta2 = $clipped - $mean;
        $m2 += $delta * $delta2;
    }
    my $sd = $seen > 1 ? sqrt($m2 / ($seen - 1)) : 0;
    return {
        n                    => $seen,
        mean                 => $mean,
        sd                   => $sd,
        clip_lower_value     => $lower_value,
        clip_upper_value     => $upper_value,
        lower_tail_clipped_n => $lower_tail_clipped_n,
        upper_tail_clipped_n => $upper_tail_clipped_n,
    };
}

sub read_numeric_values {
    my ($path, $col_i) = @_;
    open my $fh, '-|', "zcat '$path'" or die "Cannot read $path with zcat: $!\n";
    <$fh>;
    my @values;
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        my @v = split /\t/, $line, -1;
        my $x = numeric($v[$col_i]);
        push @values, $x if defined $x;
    }
    close $fh;
    return @values;
}

sub empirical_quantile {
    my ($sorted, $q) = @_;
    my $n = scalar @{$sorted};
    return undef unless $n;
    return $sorted->[0] if $n == 1 || $q <= 0;
    return $sorted->[$n - 1] if $q >= 1;
    my $pos = ($n - 1) * $q;
    my $lo = int($pos);
    my $hi = $lo < $n - 1 ? $lo + 1 : $lo;
    my $frac = $pos - $lo;
    return $sorted->[$lo] if $hi == $lo;
    return $sorted->[$lo] + (($sorted->[$hi] - $sorted->[$lo]) * $frac);
}

sub value {
    my ($vals, $idx, $col) = @_;
    return '' unless exists $idx->{$col};
    return $vals->[ $idx->{$col} ] // '';
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

sub validate_method {
    my ($method) = @_;
    my %allowed = map { $_ => 1 } qw(mean_sd mean_sd_clipped);
    die "Unsupported --method $method. Allowed values: mean_sd, mean_sd_clipped\n"
        unless defined $method && $allowed{$method};
}

sub validate_quantile {
    my ($value, $name) = @_;
    die "--$name is required\n" unless defined $value;
    my $numeric = numeric($value);
    die "--$name must be numeric\n" unless defined $numeric;
    die "--$name must be between 0 and 1 inclusive\n"
        unless $numeric >= 0 && $numeric <= 1;
    return $numeric;
}

sub usage {
    return <<"USAGE";
Usage:
  perl standardize_diff_gwas_zscore.pl [options]

Options:
  --input FILE.tsv.gz       Differential GWAS table
  --output FILE.tsv.gz      Output with ORIG_DIFF_Z, ORIG_DIFF_P, STD_DIFF_Z, STD_DIFF_P
  --manifest FILE.tsv       Run summary with mean and sample SD
  --z-col NAME              Z-score column to standardize. Default: DIFF_Z
  --method NAME             mean_sd|mean_sd_clipped . Default: mean_sd
  --clip-lower-quantile Q   Lower winsorization quantile for mean_sd_clipped. Default: 0.001
  --clip-upper-quantile Q   Upper winsorization quantile for mean_sd_clipped. Default: 0.999
  --htsbin DIR              Directory containing bgzip/tabix. Default: ../local/bin
  --index-output            Try to bgzip-index the standardized output with tabix. Default: on
  --no-index-output         Write plain gzip output without tabix indexing
USAGE
}
