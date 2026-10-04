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

my $pvar = "$prefix.pvar.zst";
my $pvar_has_rsids = -s $pvar && pvar_has_rsids($pvar);
if (-s "$prefix.pgen" && $pvar_has_rsids && -s "$prefix.psam"
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
my $psam_name = 'hg38_corrected.psam';
my $pgen_url = official_link($html, $pgen_name);
my $pvar_url = official_link($html,
    "chr${chr}_hg38_rs_noannot.pvar.zst",
    "chr${chr}_hg38_rs.pvar.zst");
my $psam_url = official_link($html, $psam_name);
download($curl, $pgen_url, "$prefix.pgen.zst") unless -s "$prefix.pgen";
if (!$pvar_has_rsids) {
    my $staged_pvar = "$prefix.rsid.pvar.zst";
    unlink $staged_pvar if -s $staged_pvar && !pvar_has_rsids($staged_pvar);
    download($curl, $pvar_url, $staged_pvar);
    die "Downloaded PVAR has no rsIDs: $staged_pvar\n"
        unless pvar_has_rsids($staged_pvar);
    if (-e $pvar) {
        my $backup = "$prefix.no_rsid.pvar.zst";
        my $suffix = 1;
        $backup = "$prefix.no_rsid.$suffix.pvar.zst" while -e $backup && $suffix++;
        rename $pvar, $backup or die "Cannot preserve $pvar as $backup: $!\n";
        print "[backup] Preserved PVAR without rsIDs: $backup\n";
        unless (rename $staged_pvar, $pvar) {
            my $error = $!;
            rename $backup, $pvar or warn "Cannot restore $pvar from $backup: $!\n";
            die "Cannot install rsID PVAR $pvar: $error\n";
        }
    } else {
        rename $staged_pvar, $pvar or die "Cannot install rsID PVAR $pvar: $!\n";
    }
}
my $shared_psam = File::Spec->catfile($output_dir, $psam_name);
download($curl, $psam_url, $shared_psam) unless -s $shared_psam;
copy($shared_psam, "$prefix.psam") or die "Cannot prepare $prefix.psam: $!\n"
    unless -s "$prefix.psam";
if (!-s "$prefix.pgen") {
    print "[prepare] Decompressing official chr$chr phased genotypes\n";
    my $tmp_pgen = "$prefix.pgen.tmp.$$";
    my ($plink2_input, $plink2_output) = ("$prefix.pgen.zst", $tmp_pgen);
    if ($^O eq 'cygwin' && $plink2 =~ /\.exe$/i) {
        # Native Windows PLINK2 cannot open Cygwin /cygdrive/... paths.
        $plink2_input = Cygwin::posix_to_win_path($plink2_input, 1);
        $plink2_output = Cygwin::posix_to_win_path($plink2_output, 1);
    }
    system($plink2, '--zst-decompress', $plink2_input, $plink2_output) == 0
        or die "PLINK2 could not decompress $prefix.pgen.zst\n";
    die "Decompressed PGEN is empty\n" unless -s $tmp_pgen;
    rename $tmp_pgen, "$prefix.pgen" or die "Cannot install $prefix.pgen: $!\n";
}
die "Official hg38 chromosome fileset is incomplete: $prefix\n"
    unless -s "$prefix.pgen" && pvar_has_rsids($pvar) && -s "$prefix.psam";
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
    my ($html, @filenames) = @_;
    for my $filename (@filenames) {
        # Match the URL filename: a displayed .pvar.zst link can point to a .log.
        my ($url) = $html =~ /href="([^"]*\/\Q$filename\E(?:\?[^"]*)?)"/i;
        next unless $url;
        $url =~ s/&amp;/&/g;
        die "Unexpected download host for $filename: $url\n"
            unless $url =~ m{^https://www\.dropbox\.com/};
        print "[source] Using official rsID PVAR $filename\n" if @filenames > 1;
        return $url;
    }
    die "Official PLINK2 resources page has no link for "
        . join(' or ', @filenames) . "\n";
}

sub pvar_has_rsids {
    my ($path) = @_;
    return 0 unless -s $path;
    my $zstdcat = $^O eq 'cygwin' ? '/usr/bin/zstdcat' : 'zstdcat';
    open my $pipe, '-|', $zstdcat, $path
        or die "Cannot inspect PVAR $path with $zstdcat: $!\n";
    my ($rows, $rsids) = (0, 0);
    while (my $line = <$pipe>) {
        next if $line =~ /^#/;
        my @fields = split /\t/, $line, 4;
        ++$rows;
        ++$rsids if defined $fields[2] && $fields[2] =~ /^rs\d+$/i;
        last if $rows >= 1000;
    }
    close $pipe; # Stopping after a sample can give the decompressor SIGPIPE.
    return $rows > 0 && $rsids >= $rows * 0.05;
}

sub download {
    my ($curl_bin, $url, $dest) = @_;
    return if -s $dest;
    my $part = "$dest.part";
    # A native Windows curl cannot open Cygwin's /cygdrive/... output path.
    my $curl_part = $^O eq 'cygwin' && $curl_bin =~ /\.exe$/i
        ? Cygwin::posix_to_win_path($part, 1)
        : $part;
    print "[download] " . basename($dest) . " from official PLINK2 resource\n";
    system($curl_bin, '-L', '--fail', '--retry', '3', '--retry-delay', '3',
        '--continue-at', '-', '-o', $curl_part, $url) == 0
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
