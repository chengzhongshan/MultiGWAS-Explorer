#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Spec;
use Cwd qw(abs_path);
use JSON::PP qw(decode_json encode_json);
use IO::Uncompress::Gunzip qw($GunzipError);
use lib File::Spec->catdir($Bin, File::Spec->updir());
use HTSToolResolver qw(resolve_hts_tool external_tool_path);

my ($input, $output, $bgzip_bin, $tabix_bin) = ('', '', '', '');
GetOptions(
    'input=s' => \$input,
    'output=s' => \$output,
    'bgzip-bin=s' => \$bgzip_bin,
    'tabix-bin=s' => \$tabix_bin,
) or die "Invalid wide-table indexing options\n";
die "--input and --output are required\n" unless length($input) && length($output);
die "Input not found: $input\n" unless -s $input;
$input = abs_path($input) || $input;
my @source_stat = stat($input);
my $manifest = "$output.manifest.json";
if (-s $output && -s "$output.tbi" && -s $manifest) {
    open my $mf, '<', $manifest or die "Cannot read $manifest: $!\n";
    local $/;
    my $old = eval { decode_json(<$mf>) };
    close $mf;
    if ($old && $old->{source} eq $input
        && $old->{source_size} == $source_stat[7]
        && $old->{source_mtime} == $source_stat[9]) {
        print "[tabix] Reusing indexed merged-wide GWAS: $output\n";
        exit 0;
    }
}

my $bgzip = resolve_hts_tool('bgzip', explicit => $bgzip_bin, start_dir => $Bin);
my $tabix = resolve_hts_tool('tabix', explicit => $tabix_bin, start_dir => $Bin);
die "Native bgzip and tabix are required for merged-wide indexing; put them in local/bin or PATH\n"
    unless $bgzip && $tabix;
make_path(dirname($output)) unless -d dirname($output);
my $stem = "$output.tmp.$$";
my $rows_file = "$stem.rows.tsv";
my $sorted_file = "$stem.sorted.tsv";
my $indexed_file = "$stem.bgz";
my $in = IO::Uncompress::Gunzip->new($input, MultiStream => 1)
    or die "Cannot read $input: $GunzipError\n";
my $header = <$in>;
die "Input is empty: $input\n" unless defined $header;
$header =~ s/[\r\n]+$//;
my @cols = split /\t/, $header, -1;
die "Merged-wide input must begin with CHR and BP columns\n"
    unless @cols >= 3 && $cols[0] eq 'CHR' && $cols[1] eq 'BP';
open my $rows_out, '>', $rows_file or die "Cannot write $rows_file: $!\n";
my ($rows, $bad_rows, $unsorted) = (0, 0, 0);
my ($last_chr, $last_bp) = ('', 0);
my %seen_chr;
while (my $line = <$in>) {
    my ($chr, $bp) = split /\t/, $line, 3;
    if (!defined($chr) || !length($chr) || !defined($bp)
        || $bp !~ /^\d+$/ || $bp < 1) {
        $bad_rows++;
        next;
    }
    if ($chr ne $last_chr) {
        $unsorted = 1 if $seen_chr{$chr};
        $seen_chr{$chr} = 1;
        $last_chr = $chr;
        $last_bp = 0;
    }
    $unsorted = 1 if $bp < $last_bp;
    $last_bp = $bp;
    print {$rows_out} $line;
    print {$rows_out} "\n" unless $line =~ /\n$/;
    $rows++;
}
close $in or die "Cannot finish reading $input: $GunzipError\n";
close $rows_out or die "Cannot close $rows_file: $!\n";
die "No coordinate-valid rows in $input\n" unless $rows;

my $data_file = $rows_file;
if ($unsorted) {
    print "[tabix] Sorting merged-wide rows by chromosome and BP for indexing\n";
    local $ENV{LC_ALL} = 'C';
    my $sort_bin = $^O =~ /cygwin/i ? '/usr/bin/sort' : 'sort';
    system($sort_bin, '-k1,1V', '-k2,2n', '-o', $sorted_file, $rows_file) == 0
        or die "Cannot sort merged-wide rows for tabix\n";
    unlink $rows_file or die "Cannot remove $rows_file: $!\n";
    $data_file = $sorted_file;
}

pipe(my $bgzip_reader, my $bgzip_out) or die "Cannot create bgzip pipe: $!\n";
my $bgzip_pid = fork();
die "Cannot fork bgzip: $!\n" unless defined $bgzip_pid;
if ($bgzip_pid == 0) {
    close $bgzip_out;
    open STDIN, '<&', $bgzip_reader or die "Cannot connect bgzip stdin: $!\n";
    open STDOUT, '>', $indexed_file or die "Cannot open $indexed_file: $!\n";
    exec { $bgzip } $bgzip, '-c';
    die "Cannot execute $bgzip: $!\n";
}
close $bgzip_reader;
print {$bgzip_out} "$header\n";
open my $rows_in, '<', $data_file or die "Cannot read $data_file: $!\n";
while (read($rows_in, my $buffer, 1024 * 1024)) {
    print {$bgzip_out} $buffer or die "Cannot stream rows to bgzip: $!\n";
}
close $rows_in or die "Cannot close $data_file: $!\n";
close $bgzip_out or die "bgzip failed for $indexed_file\n";
waitpid($bgzip_pid, 0);
die "bgzip failed for $indexed_file\n" unless $? == 0;
unlink $data_file or die "Cannot remove $data_file: $!\n";
die "bgzip produced no output: $indexed_file\n" unless -s $indexed_file;
system($tabix, '-f', '-s', '1', '-b', '2', '-e', '2', '-S', '1',
    external_tool_path($tabix, $indexed_file)) == 0
    or die "tabix failed for $indexed_file\n";
die "tabix index missing: $indexed_file.tbi\n" unless -s "$indexed_file.tbi";
for my $pair ([$indexed_file, $output], ["$indexed_file.tbi", "$output.tbi"]) {
    unlink $pair->[1] if -e $pair->[1];
    rename $pair->[0], $pair->[1] or die "Cannot install $pair->[1]: $!\n";
}
my $tmp_manifest = "$manifest.tmp.$$";
open my $mf, '>', $tmp_manifest or die "Cannot write $tmp_manifest: $!\n";
print {$mf} encode_json({
    source => $input, source_size => $source_stat[7],
    source_mtime => $source_stat[9], rows_indexed => $rows,
    rows_skipped_bad_coordinate => $bad_rows, sorted_before_indexing => ($unsorted ? 1 : 0),
}), "\n";
close $mf or die "Cannot close $tmp_manifest: $!\n";
unlink $manifest if -e $manifest;
rename $tmp_manifest, $manifest or die "Cannot install $manifest: $!\n";
print "[tabix] Indexed $rows merged-wide rows: $output; skipped bad coordinates=$bad_rows\n";
