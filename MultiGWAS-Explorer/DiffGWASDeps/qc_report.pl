#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use IO::Uncompress::Gunzip qw($GunzipError);
use JSON::PP qw(decode_json encode_json);

my ($input, $config, $outdir, $snps, $reservoir_size, $allow_filtered) = ('', '', 'qc', '', 200_000, 0);
GetOptions(
    'input=s'          => \$input,
    'config=s'         => \$config,
    'output-dir=s'     => \$outdir,
    'nominated-snps=s' => \$snps,
    'reservoir-size=i' => \$reservoir_size,
    'allow-filtered-input!' => \$allow_filtered,
) or die usage();
die usage() unless $input && $config;
die "--reservoir-size must be at least 1000\n" if $reservoir_size < 1000;
if (!$allow_filtered && $input =~ /(?:p_lt_|significant|subset)/i) {
    die "Refusing to estimate genomic inflation from a P-selected input ($input). " .
        "Use the unfiltered genome-wide wide table, or pass --allow-filtered-input only for debugging; " .
        "filtered lambda values are not publication-valid.\n";
}

open my $cfh, '<', $config or die "Cannot read $config: $!\n";
my $cfg_text = do { local $/; <$cfh> };
my $cfg = decode_json($cfg_text);
close $cfh;
make_path($outdir) unless -d $outdir;
die "Cannot create output directory $outdir\n" unless -d $outdir;

my @prefixes = map { $_->{prefix} } @{ $cfg->{pairs} || [] };
die "Config has no pair prefixes\n" unless @prefixes;
my %wanted = map { $_ => 1 } grep { length } split /\s*,\s*/, $snps;

my $fh = open_reader($input);
my $header = <$fh> // die "Empty input: $input\n";
chomp $header;
$header =~ s/\r$//;
my @cols = split /\t/, $header, -1;
my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
for my $required (qw(CHR BP SNP)) {
    die "Required column $required is absent from $input\n" unless exists $idx{$required};
}

my @tracks;
for my $prefix (@prefixes) {
    for my $group (1, 2) {
        my $id = "${prefix}_GROUP$group";
        my ($b, $s, $p) = map { "${id}_$_" } qw(BETA SE P);
        push @tracks, { id => $id, beta => $b, se => $s, p => $p }
          if exists $idx{$b} && exists $idx{$s};
    }
    my $z = "${prefix}_STD_DIFF_Z";
    $z = "${prefix}_DIFF_Z" unless exists $idx{$z};
    push @tracks, { id => "${prefix}_DIFFERENTIAL", z => $z,
                    p => (exists $idx{"${prefix}_STD_DIFF_P"} ? "${prefix}_STD_DIFF_P" : "${prefix}_DIFF_P") }
      if exists $idx{$z};
}
die "No association or differential tracks were detected\n" unless @tracks;

my %state = map { $_->{id} => { n => 0, sample => [] } } @tracks;
my @hits;
my ($rows, $autosomal, $chrx) = (0, 0, 0);
srand(20260819);
while (my $line = <$fh>) {
    chomp $line;
    $line =~ s/\r$//;
    next unless length $line;
    my @v = split /\t/, $line, -1;
    $rows++;
    my $chr = $v[$idx{CHR}] // '';
    $chr =~ s/^chr//i;
    $chr =~ /^(?:X|23)$/i ? $chrx++ : $autosomal++ if $chr =~ /^(?:[1-9]|1\d|2[0-3]|X)$/i;

    for my $t (@tracks) {
        my $z;
        if ($t->{z}) {
            $z = number($v[$idx{$t->{z}}]);
        } else {
            my $b = number($v[$idx{$t->{beta}}]);
            my $s = number($v[$idx{$t->{se}}]);
            $z = $b / $s if defined $b && defined $s && $s > 0;
        }
        next unless defined $z;
        reservoir_add($state{$t->{id}}, $z * $z, $reservoir_size);
    }

    if (%wanted && $wanted{$v[$idx{SNP}] // ''}) {
        my %row = map { ($_, ($v[$idx{$_}] // '')) } grep {
            /^(?:CHR|BP|A1|A2|SNP)$/ || /_(?:GROUP[12]_(?:BETA|SE|P|FRQ_A|FRQ_U|INFO)|DIFF_(?:BETA|SE|P)|STD_DIFF_(?:Z|P))$/
        } @cols;
        push @hits, \%row;
    }
}
close $fh;

open my $lfh, '>', "$outdir/genomic_inflation.tsv" or die $!;
print {$lfh} "TRACK\tN_VALID\tN_RESERVOIR\tMEDIAN_CHISQ\tLAMBDA_GC\tMETHOD\n";
my @lambda;
for my $t (@tracks) {
    my $s = $state{$t->{id}};
    my @x = sort { $a <=> $b } @{ $s->{sample} };
    my $median = median(\@x);
    my $lambda = defined $median ? $median / 0.454936423119573 : undef;
    print {$lfh} join("\t", $t->{id}, $s->{n}, scalar(@x), fmt($median), fmt($lambda),
        $s->{n} > @x ? 'deterministic_seed_reservoir_estimate' : 'exact'), "\n";
    push @lambda, { track => $t->{id}, n_valid => $s->{n}, n_reservoir => scalar(@x),
                   median_chisq => $median, lambda_gc => $lambda };
}
close $lfh;

my @hit_cols = grep { /^(?:CHR|BP|A1|A2|SNP)$/ || /_(?:GROUP[12]_(?:BETA|SE|P|FRQ_A|FRQ_U|INFO)|DIFF_(?:BETA|SE|P)|STD_DIFF_(?:Z|P))$/ } @cols;
open my $hfh, '>', "$outdir/nominated_variants.tsv" or die $!;
print {$hfh} join("\t", @hit_cols, qw(QC_MAF_MIN QC_INFO_MIN QC_FLAGS)), "\n";
for my $r (@hits) {
    my (@mafs, @infos, @flags);
    for my $c (@hit_cols) {
        my $x = number($r->{$c});
        if ($c =~ /_FRQ_[AU]$/ && defined $x && $x >= 0 && $x <= 1) { push @mafs, $x > .5 ? 1-$x : $x }
        if ($c =~ /_INFO$/ && defined $x) { push @infos, $x }
    }
    my $maf = @mafs ? min(@mafs) : undef;
    my $info = @infos ? min(@infos) : undef;
    push @flags, 'MAF_AT_OR_BELOW_0.01' if defined $maf && $maf <= .01;
    push @flags, 'MAF_MISSING' unless defined $maf;
    push @flags, 'INFO_BELOW_0.8' if defined $info && $info < .8;
    push @flags, 'INFO_MISSING' unless defined $info;
    print {$hfh} join("\t", (map { $r->{$_} // '' } @hit_cols), fmt($maf), fmt($info), (@flags ? join(';', @flags) : 'PASS')), "\n";
}
close $hfh;

my $summary = {
    schema_version => 1, input => $input, config => $config, rows_scanned => $rows,
    coordinate_rows_autosomal => $autosomal, coordinate_rows_chrX => $chrx,
    nominated_snps_requested => [sort keys %wanted], nominated_rows_found => scalar(@hits),
    genomic_inflation => \@lambda,
    cautions => [
        'Lambda GC is descriptive and does not replace LD-score-regression intercepts.',
        'A summary-statistic difference test is not a substitute for an individual-level genotype-by-sex interaction model.',
        'LD and conditional independence require an ancestry-matched reference panel or individual-level genotypes and are not inferred from physical distance.',
    ],
};
open my $jfh, '>', "$outdir/qc_summary.json" or die $!;
print {$jfh} JSON::PP->new->canonical->pretty->encode($summary);
close $jfh;
warn "QC report written to $outdir ($rows rows scanned)\n";

sub open_reader {
    my ($path) = @_;
    if ($path =~ /\.gz$/i) {
        my $x = IO::Uncompress::Gunzip->new($path)
          or die "Cannot open $path: $GunzipError\n";
        return $x;
    }
    open my $x, '<', $path or die "Cannot open $path: $!\n";
    return $x;
}
sub number {
    my ($x) = @_;
    return undef unless defined $x && $x =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?$/;
    return 0 + $x;
}
sub reservoir_add {
    my ($s, $x, $cap) = @_;
    $s->{n}++;
    if (@{$s->{sample}} < $cap) { push @{$s->{sample}}, $x; return }
    my $j = int(rand($s->{n}));
    $s->{sample}[$j] = $x if $j < $cap;
}
sub median {
    my ($x) = @_;
    return undef unless @$x;
    my $n = @$x;
    return $n % 2 ? $x->[int($n/2)] : ($x->[$n/2-1] + $x->[$n/2]) / 2;
}
sub min {
    my $m = $_[0];
    for my $x (@_[1..$#_]) { $m = $x if $x < $m }
    return $m;
}
sub fmt { defined $_[0] ? sprintf('%.10g', $_[0]) : '' }
sub usage { "Usage: perl qc_report.pl --input UNFILTERED_WIDE.tsv[.gz] --config SPEC.json --output-dir DIR [--nominated-snps rs1,rs2] [--reservoir-size 200000] [--allow-filtered-input]\n" }
