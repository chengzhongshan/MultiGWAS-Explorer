#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(basename dirname);
use File::Copy qw(copy);
use File::Path qw(make_path);
use File::Spec;
use IO::Uncompress::Unzip qw($UnzipError);

my ($chr, $output_dir, $plink2, $curl) = ('', '', '', '');
GetOptions(
    'chr=s' => \$chr,
    'output-dir=s' => \$output_dir,
    'plink2=s' => \$plink2,
    'curl=s' => \$curl,
) or die "Invalid 1000 Genomes reference options\n";
die "--chr and --output-dir are required\n" unless length($chr) && length($output_dir);
$chr =~ s/^chr//i;
$chr = 'X' if $chr eq '23';
die "--chr must be 1-22 or X\n"
    unless $chr eq 'X' || ($chr =~ /^\d+$/ && $chr >= 1 && $chr <= 22);
$curl ||= $^O =~ /cygwin|MSWin32/i ? 'curl.exe' : 'curl';
make_path($output_dir) unless -d $output_dir;
my $prefix = File::Spec->catfile($output_dir, "chr${chr}_hg38");
$plink2 ||= File::Spec->catfile(dirname($output_dir), 'plink2_bin',
    $^O =~ /cygwin|MSWin32/i ? 'plink2.exe' : 'plink2');
chmod 0755, $plink2 if -f $plink2 && !-x $plink2;

if (-s "$prefix.pgen" && -s "$prefix.pvar.zst" && -s "$prefix.psam"
    && -f $plink2 && -x $plink2) {
    print "[reuse] Official PLINK2 1000 Genomes Phase 3 hg38 chr$chr: $prefix\n";
    exit 0;
}

my $resources = 'https://www.cog-genomics.org/plink/2.0/resources';
my $html = fetch_page($curl, $resources);
if (!(-f $plink2 && -x $plink2)) {
    my $downloads = fetch_page($curl, 'https://www.cog-genomics.org/plink/2.0/');
    my $zip_pattern = $^O =~ /cygwin|MSWin32/i
        ? qr/plink2_win64_\d+\.zip/i
        : qr/plink2_linux_x86_64_\d+\.zip/i;
    my ($zip_url) = $downloads =~ /href="([^"]*$zip_pattern[^"]*)"/i;
    die "Cannot find a native PLINK2 binary on the official download page\n"
        unless $zip_url;
    $zip_url =~ s/&amp;/&/g;
    make_path(dirname($plink2)) unless -d dirname($plink2);
    my $zip_path = "$plink2.download.zip";
    download($curl, $zip_url, $zip_path);
    extract_plink2($zip_path, $plink2);
    unlink $zip_path;
}
system($plink2, '--version') == 0
    or die "PLINK2 executable does not run: $plink2\n";

my $pgen_name = "chr${chr}_hg38.pgen.zst";
my $pvar_name = "chr${chr}_hg38_rs_noannot.pvar.zst";
my $psam_name = 'hg38_corrected.psam';
my $pgen_url = official_link($html, $pgen_name);
my $pvar_url = official_link($html, $pvar_name);
my $psam_url = official_link($html, $psam_name);
download($curl, $pgen_url, "$prefix.pgen.zst") unless -s "$prefix.pgen";
download($curl, $pvar_url, "$prefix.pvar.zst") unless -s "$prefix.pvar.zst";
my $shared_psam = File::Spec->catfile($output_dir, $psam_name);
download($curl, $psam_url, $shared_psam) unless -s $shared_psam;
copy($shared_psam, "$prefix.psam") or die "Cannot prepare $prefix.psam: $!\n"
    unless -s "$prefix.psam";
if (!-s "$prefix.pgen") {
    print "[prepare] Decompressing official chr$chr phased genotypes\n";
    my $tmp_pgen = "$prefix.pgen.tmp.$$";
    system($plink2, '--zst-decompress', "$prefix.pgen.zst", $tmp_pgen) == 0
        or die "PLINK2 could not decompress $prefix.pgen.zst\n";
    die "Decompressed PGEN is empty\n" unless -s $tmp_pgen;
    rename $tmp_pgen, "$prefix.pgen" or die "Cannot install $prefix.pgen: $!\n";
}
die "Official hg38 chromosome fileset is incomplete: $prefix\n"
    unless -s "$prefix.pgen" && -s "$prefix.pvar.zst" && -s "$prefix.psam";
print "[ready] PLINK2 1000 Genomes Phase 3 GRCh38/hg38 chr$chr: $prefix\n";

sub fetch_page {
    my ($curl_bin, $url) = @_;
    open my $pipe, '-|', $curl_bin, '-LfsS', $url
        or die "Cannot fetch $url: $!\n";
    local $/;
    my $body = <$pipe> // '';
    close $pipe or die "Cannot fetch $url\n";
    die "Official download page is empty: $url\n" unless length $body;
    return $body;
}

sub official_link {
    my ($html, $filename) = @_;
    my ($url) = $html =~ /href="([^"]*\Q$filename\E[^"]*)"/i;
    die "Official PLINK2 resources page has no link for $filename\n" unless $url;
    $url =~ s/&amp;/&/g;
    die "Unexpected download host for $filename: $url\n"
        unless $url =~ m{^https://www\.dropbox\.com/};
    return $url;
}

sub download {
    my ($curl_bin, $url, $dest) = @_;
    return if -s $dest;
    my $part = "$dest.part";
    print "[download] " . basename($dest) . " from official PLINK2 resource\n";
    system($curl_bin, '-L', '--fail', '--retry', '3', '--retry-delay', '3',
        '--continue-at', '-', '-o', $part, $url) == 0
        or die "Download failed for $dest\n";
    die "Downloaded file is empty: $dest\n" unless -s $part;
    rename $part, $dest or die "Cannot install $dest: $!\n";
}

sub extract_plink2 {
    my ($zip_path, $dest) = @_;
    my $zip = IO::Uncompress::Unzip->new($zip_path)
        or die "Cannot open PLINK2 ZIP $zip_path: $UnzipError\n";
    my $expected = $^O =~ /cygwin|MSWin32/i ? 'plink2.exe' : 'plink2';
    my $found = 0;
    for (;;) {
        my $name = $zip->getHeaderInfo->{Name} // '';
        if (lc(basename($name)) eq lc($expected)) {
            my $tmp = "$dest.tmp.$$";
            open my $out, '>', $tmp or die "Cannot write $tmp: $!\n";
            while (my $n = $zip->read(my $buffer, 1024 * 1024)) {
                print {$out} $buffer or die "Cannot write $tmp: $!\n";
            }
            close $out or die "Cannot close $tmp: $!\n";
            chmod 0755, $tmp or die "Cannot make $tmp executable: $!\n";
            rename $tmp, $dest or die "Cannot install $dest: $!\n";
            $found = 1;
            last;
        }
        last unless $zip->nextStream() == 1;
    }
    die "PLINK2 binary was not found in $zip_path\n" unless $found && -s $dest;
}
