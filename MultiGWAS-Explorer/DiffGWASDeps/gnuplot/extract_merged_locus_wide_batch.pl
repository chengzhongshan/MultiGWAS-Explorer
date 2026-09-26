#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Spec;
use JSON::PP qw(decode_json);
use IO::Compress::Gzip qw($GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);

my ($input, $indexed_input, $tabix_bin, $output_dir, $window_bp) = ('', '', '', '', '');
my @target_args;
GetOptions(
    'input=s'      => \$input,
    'indexed-input=s' => \$indexed_input,
    'tabix-bin=s'  => \$tabix_bin,
    'output-dir=s' => \$output_dir,
    'window-bp=s'  => \$window_bp,
    'target=s@'    => \@target_args,
) or die "Invalid locus extraction options\n";
die "--input, --output-dir, --window-bp, and at least one --target are required\n"
    unless length($input) && length($output_dir) && length($window_bp) && @target_args;
die "Input not found: $input\n" unless -s $input;
my @source_stat = stat($input);
die "Indexed input or tabix index is missing: $indexed_input\n"
    if $indexed_input && !(-s $indexed_input && -s "$indexed_input.tbi");
die "--window-bp must be positive\n"
    unless $window_bp =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i && $window_bp > 0;
make_path($output_dir) unless -d $output_dir;

my (%targets_by_chr, @targets, %seen);
for my $arg (@target_args) {
    my $spec = eval { decode_json($arg) };
    die "Invalid --target JSON: $arg\n" unless ref($spec) eq 'HASH';
    my ($snp, $chr, $bp) = @{$spec}{qw(snp chr bp)};
    die "Invalid --target coordinate: $arg\n"
        unless defined($snp) && length($snp) && defined($chr)
            && defined($bp) && $bp =~ /^\d+$/;
    $chr = normalize_chr($chr);
    die "Invalid target coordinate: $arg\n" unless length($chr) && $bp > 0;
    next if $seen{uc($snp)}++;
    my $stem = 'gnuplot_locus_' . safe_name($snp)
        . '_window_' . safe_name($window_bp) . '.wide';
    my $target = {
        snp => $snp, chr => $chr, bp => 0 + $bp,
        start => $bp - $window_bp > 1 ? $bp - $window_bp : 1,
        end => $bp + $window_bp,
        data => File::Spec->catfile($output_dir, "$stem.tsv.gz"),
        manifest => File::Spec->catfile($output_dir, "$stem.manifest.tsv"),
        rows_written => 0, target_found => 0,
    };
    push @targets, $target;
    push @{ $targets_by_chr{$chr} }, $target;
}
die "No distinct target SNPs were supplied\n" unless @targets;

my $in = IO::Uncompress::Gunzip->new($input, MultiStream => 1)
    or die "Cannot read $input: $GunzipError\n";
my $header = <$in>;
die "Input is empty: $input\n" unless defined $header;
my $header_text = $header;
$header_text =~ s/[\r\n]+$//;
my @cols = split /\t/, $header_text, -1;
my %idx = map { $cols[$_] => $_ } 0 .. $#cols;
for my $required (qw(CHR BP SNP)) {
    die "Required column $required is missing from $input\n" unless exists $idx{$required};
}
for my $target (@targets) {
    $target->{tmp_data} = "$target->{data}.tmp.$$";
    my $out = IO::Compress::Gzip->new($target->{tmp_data})
        or die "Cannot write $target->{tmp_data}: $GzipError\n";
    print {$out} $header_text, "\n";
    $target->{out} = $out;
}

my $rows_read = 0;
my $record = sub {
    my ($line, $target_list) = @_;
    $rows_read++;
    my @f = split /\t/, $line, -1;
    my $chr = normalize_chr($f[$idx{CHR}] // '');
    my $bp = $f[$idx{BP}] // '';
    return unless $target_list && $bp =~ /^\d+$/ && $bp > 0;
    my $snp = $f[$idx{SNP}] // '';
    $snp =~ s/[\r\n]+$//;
    for my $target (@$target_list) {
        next if $bp < $target->{start} || $bp > $target->{end};
        print { $target->{out} } $line;
        print { $target->{out} } "\n" unless $line =~ /\n$/;
        $target->{rows_written}++;
        $target->{target_found} = 1 if lc($snp) eq lc($target->{snp});
    }
};
if ($indexed_input) {
    close $in or die "Cannot close $input: $GunzipError\n";
    my $tabix = find_tool('tabix', $tabix_bin || $ENV{TABIX_BIN});
    die "Native tabix is required to query $indexed_input\n" unless $tabix;
    open my $contigs, '-|', $tabix, '-l', $indexed_input
        or die "Cannot list tabix contigs in $indexed_input: $!\n";
    my %indexed_chr;
    while (my $contig = <$contigs>) {
        chomp $contig;
        my $normalized = normalize_chr($contig);
        push @{ $indexed_chr{$normalized} }, $contig if length $normalized;
    }
    close $contigs or die "Cannot list tabix contigs in $indexed_input\n";
    for my $target (@targets) {
        for my $contig (@{ $indexed_chr{$target->{chr}} || [] }) {
            my $region = "$contig:$target->{start}-$target->{end}";
            open my $query, '-|', $tabix, $indexed_input, $region
                or die "Cannot query $region in $indexed_input: $!\n";
            while (my $line = <$query>) {
                $record->($line, [$target]);
            }
            close $query or die "tabix query failed: $region\n";
        }
    }
}
else {
while (my $line = <$in>) {
    my ($chr_text) = split /\t/, $line, 2;
    my $chr = normalize_chr($chr_text);
    $record->($line, $targets_by_chr{$chr});
}
close $in or die "Cannot finish reading $input: $GunzipError\n";
}
for my $target (@targets) {
    close $target->{out} or die "Cannot finish $target->{tmp_data}: $GzipError\n";
    die "Target SNP $target->{snp} was not found at $target->{chr}:$target->{bp} in $input\n"
        unless $target->{target_found};
    replace_file($target->{tmp_data}, $target->{data});
    my $tmp_manifest = "$target->{manifest}.tmp.$$";
    open my $mf, '>', $tmp_manifest or die "Cannot write $tmp_manifest: $!\n";
    print {$mf} "METRIC\tVALUE\n";
    for my $entry (
        ['target_snp', $target->{snp}],
        ['target_chr', $target->{chr}],
        ['target_bp', $target->{bp}],
        ['window_bp', $window_bp],
        ['rows_written', $target->{rows_written}],
        ['source', $input],
        ['source_size', $source_stat[7]],
        ['source_mtime', $source_stat[9]],
        ['access_mode', ($indexed_input ? 'TABIX' : 'STREAM')],
    ) {
        print {$mf} join("\t", @$entry), "\n";
    }
    close $mf or die "Cannot close $tmp_manifest: $!\n";
    replace_file($tmp_manifest, $target->{manifest});
    print "[locus] $target->{snp}: $target->{rows_written} rows -> $target->{data}\n";
}
print $indexed_input
    ? "[locus] Queried " . scalar(@targets) . " indexed windows ($rows_read rows returned)\n"
    : "[locus] Scanned $rows_read source rows once for " . scalar(@targets) . " target windows\n";

sub normalize_chr {
    my ($chr) = @_;
    $chr //= '';
    $chr =~ s/^chr//i;
    $chr = 23 if uc($chr) eq 'X';
    $chr = 24 if uc($chr) eq 'Y';
    return $chr =~ /^\d+$/ && $chr >= 1 && $chr <= 24 ? '' . (0 + $chr) : '';
}

sub safe_name {
    my ($text) = @_;
    $text =~ s/[^A-Za-z0-9._-]+/_/g;
    $text =~ s/^_+|_+$//g;
    return length($text) ? $text : 'item';
}

sub replace_file {
    my ($from, $to) = @_;
    unlink $to if -e $to;
    rename $from, $to or die "Cannot replace $to: $!\n";
}

sub find_tool {
    my ($name, $explicit) = @_;
    if ($explicit) {
        return $explicit if -f $explicit && -x $explicit;
        die "Configured $name is unavailable: $explicit\n";
    }
    my @suffixes = $^O =~ /^(?:cygwin|MSWin32)$/i ? ('.exe', '') : ('');
    my @dirs = ($Bin);
    my $ancestor = $Bin;
    for (1 .. 5) {
        $ancestor = dirname($ancestor);
        push @dirs, File::Spec->catdir($ancestor, 'local', 'bin');
    }
    push @dirs, '/usr/bin' if $^O =~ /cygwin/i;
    push @dirs, File::Spec->path();
    my %seen;
    for my $dir (grep { defined($_) && length($_) && !$seen{$_}++ } @dirs) {
        for my $suffix (@suffixes) {
            my $candidate = File::Spec->catfile($dir, $name . $suffix);
            return $candidate if -f $candidate && -x $candidate;
        }
    }
    return;
}
