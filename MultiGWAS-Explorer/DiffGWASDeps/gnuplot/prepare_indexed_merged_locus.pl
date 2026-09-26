#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);
use Cwd qw(abs_path);
use Digest::SHA qw(sha1_hex);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use File::Spec;
use JSON::PP qw(encode_json);
use IO::Uncompress::Gunzip qw($GunzipError);

my ($input, $target_snp, $window_bp, $output_dir, $cache_dir, $target_chr, $target_bp) = ('') x 7;
GetOptions(
    'input=s'       => \$input,
    'target-snp=s'  => \$target_snp,
    'window-bp=s'   => \$window_bp,
    'output-dir=s'  => \$output_dir,
    'cache-dir=s'   => \$cache_dir,
    'target-chr=s'  => \$target_chr,
    'target-bp=i'   => \$target_bp,
) or die "Invalid indexed merged-locus options\n";
die "--input, --target-snp, and --window-bp are required\n"
    unless length($input) && length($target_snp) && length($window_bp);
die "Merged-wide input does not exist: $input\n" unless -s $input;
die "--window-bp must be positive\n"
    unless $window_bp =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i && $window_bp > 0;
$input = abs_path($input) || $input;
$output_dir ||= dirname($input);
$cache_dir ||= File::Spec->catdir($Bin, '..', '..', 'cache', 'gnuplot_wide_index');
make_path($output_dir) unless -d $output_dir;
make_path($cache_dir) unless -d $cache_dir;

my $safe_snp = safe_name($target_snp);
my $safe_window = safe_name($window_bp);
my $data = File::Spec->catfile($output_dir,
    "gunplot_locus_${safe_snp}_window_${safe_window}.wide.tsv.gz");
my $manifest = File::Spec->catfile($output_dir,
    "gunplot_locus_${safe_snp}_window_${safe_window}.wide.manifest.tsv");
my @source_stat = stat($input);
my %old = -s $manifest ? read_manifest($manifest) : ();
my $source_matches = defined($old{source})
    && (abs_path($old{source}) || $old{source}) eq $input
    && ($old{source_size} // '') eq $source_stat[7]
    && ($old{source_mtime} // '') eq $source_stat[9];
if ($source_matches && ($old{target_snp} // '') eq $target_snp
    && ($old{target_bp} // '') =~ /^\d+$/) {
    $target_chr ||= $old{target_chr};
    $target_bp ||= $old{target_bp};
}

if (!$target_chr || !$target_bp) {
    my $fh = IO::Uncompress::Gunzip->new($input, MultiStream => 1)
        or die "Cannot read $input: $GunzipError\n";
    my $header = <$fh> // die "Empty merged-wide input: $input\n";
    $header =~ s/[\r\n]+$//;
    my @cols = split /\t/, $header, -1;
    my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
    die "Merged-wide input needs CHR, BP, and SNP columns\n"
        unless exists($idx{CHR}) && exists($idx{BP}) && exists($idx{SNP});
    while (my $line = <$fh>) {
        my @f = split /\t/, $line, -1;
        next unless lc($f[$idx{SNP}] // '') eq lc($target_snp);
        $target_chr = $f[$idx{CHR}];
        $target_bp = $f[$idx{BP}];
        last;
    }
    close $fh or die "Cannot finish reading $input: $GunzipError\n";
    die "Target SNP $target_snp was not found in $input\n"
        unless $target_chr && $target_bp;
    print "[locus] Located $target_snp at $target_chr:$target_bp in merged-wide GWAS\n";
}
die "Invalid target BP for $target_snp: $target_bp\n"
    unless $target_bp =~ /^\d+$/ && $target_bp > 0;
$target_chr =~ s/^chr//i;
$target_chr = 23 if uc($target_chr) eq 'X';
$target_chr = 24 if uc($target_chr) eq 'Y';
die "Invalid target chromosome for $target_snp: $target_chr\n"
    unless $target_chr =~ /^\d+$/ && $target_chr >= 1 && $target_chr <= 24;
$target_chr = 0 + $target_chr;

if (-s $data && $source_matches && ($old{target_snp} // '') eq $target_snp
    && ($old{window_bp} // '') == $window_bp
    && ($old{access_mode} // '') eq 'TABIX') {
    print "[locus] Reusing tabix-extracted target window: $data\n";
}
else {
    my $identity = $input;
    if ($^O =~ /cygwin/i) {
        if (open my $pipe, '-|', 'cygpath', '-m', $input) {
            my $native = <$pipe> // '';
            close $pipe;
            $native =~ s/[\r\n]+$//;
            $identity = $native if length $native;
        }
    }
    my $index_name = safe_name(basename($input)) . '.'
        . substr(sha1_hex($identity), 0, 12) . '.bgz';
    my $indexed = File::Spec->catfile($cache_dir, $index_name);
    run($^X, File::Spec->catfile($Bin, 'index_merged_wide_tabix.pl'),
        '--input', $input, '--output', $indexed);
    run($^X, File::Spec->catfile($Bin, 'extract_merged_locus_wide_batch.pl'),
        '--input', $input, '--indexed-input', $indexed,
        '--output-dir', $output_dir, '--window-bp', $window_bp,
        '--target', encode_json({
            snp => $target_snp, chr => $target_chr, bp => 0 + $target_bp,
        }));
}
die "Tabix extraction did not produce $data\n" unless -s $data && -s $manifest;
print "OUTPUT\t$data\nMANIFEST\t$manifest\nTARGET_CHR\t$target_chr\nTARGET_BP\t$target_bp\nACCESS_MODE\tTABIX\n";

sub run {
    my @cmd = @_;
    system(@cmd) == 0 or die "Command failed: @cmd\n";
}

sub read_manifest {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    my %values;
    while (my $line = <$fh>) {
        chomp $line;
        my ($key, $value) = split /\t/, $line, 2;
        $values{$key} = $value if defined($key) && defined($value);
    }
    close $fh;
    return %values;
}

sub safe_name {
    my ($value) = @_;
    $value =~ s/[^A-Za-z0-9._-]+/_/g;
    $value =~ s/^_+|_+$//g;
    return length($value) ? $value : 'item';
}
