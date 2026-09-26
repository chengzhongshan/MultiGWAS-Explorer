#!/usr/bin/env perl
use strict;
use warnings;

use File::Temp qw(tempdir);
use File::Spec;
use File::Path qw(make_path);
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);
use Text::ParseWords qw(shellwords);

my ($query, $pfile, $bfile, $plink2, $chr, $from_bp, $to_bp, $keep, $output, $populations);
my $min_r2 = 0;
my $window_kb = 1000;
my $phased = 1;
my $quiet = 0;
my $reference_build = 'GRCh37_hg19';

GetOptions(
    'query-snp=s' => \$query,
    'pfile=s' => \$pfile,
    'bfile=s' => \$bfile,
    'plink2=s' => \$plink2,
    'chr=s' => \$chr,
    'from-bp=i' => \$from_bp,
    'to-bp=i' => \$to_bp,
    'window-kb=f' => \$window_kb,
    'min-r2=f' => \$min_r2,
    'keep=s' => \$keep,
    'populations=s' => \$populations,
    'unphased' => sub { $phased = 0 },
    'quiet!' => \$quiet,
    'reference-build=s' => \$reference_build,
    'output=s' => \$output,
) or die usage();

die usage() unless defined($query) && (defined($pfile) xor defined($bfile));
die "--min-r2 must be between 0 and 1\n" unless $min_r2 >= 0 && $min_r2 <= 1;
die "--window-kb must be positive\n" unless $window_kb > 0;
die "--reference-build must be GRCh37_hg19 or GRCh38_hg38\n"
    unless $reference_build =~ /^(?:GRCh37_hg19|GRCh38_hg38)$/;
$plink2 ||= $ENV{PLINK2} || 'plink2';

my $tmp_root = defined($output) && length($output)
    ? dirname($output) : File::Spec->tmpdir();
make_path($tmp_root) unless -d $tmp_root;
# Cygwin's temporary directory can resolve inside an Anaconda installation
# that native Windows PLINK2 cannot open. Keep intermediate files beside the
# requested cache, where both programs already have access.
my $tmp = tempdir('plink2_local_ld.XXXXXX', DIR => $tmp_root, CLEANUP => 1);
my $prefix = "$tmp/ld";
my $input_prefix = defined($pfile) ? $pfile : $bfile;
my $pfile_arg = native_path($input_prefix);
my $tmp_prefix_arg = native_path($prefix);
if (defined($populations) && length($populations) && !defined($keep)) {
    my %wanted;
    for my $pop (split /,/, $populations) {
        $pop =~ s/^\s+|\s+$//g;
        $wanted{uc($pop)} = 1 if length $pop;
    }
    my $psam = -e "$input_prefix.psam" ? "$input_prefix.psam" : undef;
    if (!defined($psam) && $input_prefix =~ /_biallelic\z/) {
        my $candidate = $input_prefix . '';
        $candidate =~ s/_biallelic\z//;
        $candidate .= '.psam';
        $psam = $candidate if -e $candidate;
    }
    $psam ||= (-e "$input_prefix.fam" ? "$input_prefix.fam" : undef);
    $psam = undef unless defined($psam) && -e $psam;
    die "--populations requires a matching .psam file (use --keep for BED files)\n" unless $psam;
    $keep = "$tmp/populations.keep";
    open my $pfh, '<', $psam or die "Cannot read $psam: $!\n";
    open my $kfh, '>', $keep or die "Cannot write $keep: $!\n";
    while (my $line = <$pfh>) {
        next if $line =~ /^#/;
        chomp $line;
        $line =~ s/\r\z//;
        my @f = split /\s+/, $line;
        my ($iid, $super) = @f[0,4];
        print {$kfh} "$iid\n" if defined($super) && $wanted{uc($super)};
    }
    close $kfh; close $pfh;
    die "No samples matched --populations=$populations\n" unless -s $keep;
}
my @cmd;
if (-f $plink2) {
    # A literal executable path must not be passed through shellwords(): on
    # Cygwin, that strips backslashes from native Windows paths such as
    # C:\\tools\\plink2.exe.  Preserve existing paths verbatim; shellwords is
    # retained only for an explicit command string with arguments.
    @cmd = ($plink2);
} else {
    @cmd = shellwords($plink2);
}
if (defined $pfile) {
    push @cmd, '--pfile', $pfile_arg, 'vzs' if -e "$pfile.pvar.zst";
    push @cmd, '--pfile', $pfile_arg unless -e "$pfile.pvar.zst";
} else {
    push @cmd, '--bfile', $pfile_arg;
}
push @cmd, '--chr', $chr if defined($chr) && length($chr);
push @cmd, '--from-bp', $from_bp if defined $from_bp;
push @cmd, '--to-bp', $to_bp if defined $to_bp;
push @cmd, '--keep', native_path($keep) if defined($keep) && length($keep);
push @cmd, '--ld-snp', $query;
push @cmd, ($phased ? '--r2-phased' : '--r2-unphased'), 'allow-ambiguous-allele';
push @cmd, '--ld-window-kb', $window_kb;
push @cmd, '--ld-window-r2', $min_r2;
push @cmd, '--out', $tmp_prefix_arg;

system(@cmd);
die "PLINK2 LD calculation failed (exit $?): @cmd\n" if $? != 0;
my ($report) = grep { -f $_ } ("$prefix.vcor", "$prefix.vcor.zst");
die "PLINK2 did not produce a vcor report\n" unless $report;
if ($report =~ /\.zst\z/) {
    my $zcat = $ENV{ZSTD} || 'zstd';
    open my $in, '-|', $zcat, '-dc', $report or die "Cannot decompress $report: $!\n";
    parse_report($in, $output, $query, $min_r2, $quiet);
    close $in;
} else {
    open my $in, '<', $report or die "Cannot read $report: $!\n";
    parse_report($in, $output, $query, $min_r2, $quiet);
    close $in;
}

sub parse_report {
    my ($fh, $out_path, $ref, $threshold, $quiet_output) = @_;
    my ($header, %idx, %best);
    while (my $line = <$fh>) {
        chomp $line;
        next unless length $line;
        if (!$header) {
            $header = 1;
            my @h = split /\t/, $line, -1;
            s/\r\z// for @h;
            %idx = map { lc($h[$_]) =~ s/^#//r => $_ } 0 .. $#h;
            next;
        }
        my @f = split /\t/, $line, -1;
        s/\r\z// for @f;
        my $a = value(\@f, \%idx, qw(id_a variant_id_a));
        my $b = value(\@f, \%idx, qw(id_b variant_id_b));
        my $r2 = value(\@f, \%idx, qw(r2 phased_r2 unphased_r2));
        $a = $f[2] if !defined($a) && @f >= 7;
        $b = $f[5] if !defined($b) && @f >= 7;
        $r2 = $f[6] if !defined($r2) && @f >= 7;
        next unless defined($a) && defined($b) && defined($r2);
        # PLINK can store multiple aliases in one semicolon-delimited ID.
        # Expand them so a GWAS rsID can match any alias, and discard the
        # missing-ID sentinel instead of emitting an unusable proxy named '.'.
        my @a_ids = grep { length && $_ ne '.' } split /;/, $a;
        my @b_ids = grep { length && $_ ne '.' } split /;/, $b;
        next unless @a_ids && @b_ids;
        next unless $r2 =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i;
        next unless $r2 >= $threshold;
        my $a_is_ref = grep { lc($_) eq lc($ref) } @a_ids;
        my $b_is_ref = grep { lc($_) eq lc($ref) } @b_ids;
        my @proxies = $a_is_ref ? @b_ids : ($b_is_ref ? @a_ids : ());
        for my $proxy (@proxies) {
            next if lc($proxy) eq lc($ref);
            $best{$proxy} = 0 + $r2
                if !exists($best{$proxy}) || $r2 > $best{$proxy};
        }
    }
    if (defined $out_path) {
        make_path(dirname($out_path)) unless -d dirname($out_path);
    }
    open my $out, '>', $out_path or die "Cannot write $out_path: $!\n" if defined $out_path;
    print {$out} join("\t", qw(query_snp proxy_snp ld_population proxy_r2 source reference_panel reference_build ld_method)), "\n" if $out;
    if (!$quiet_output) {
        print "LD_SNPS\t", join(',', sort { $best{$b} <=> $best{$a} || $a cmp $b } keys %best), "\n";
        print "LD_R2_PAIRS\t", join(',', map { $_ . ':' . sprintf('%.6g', $best{$_}) } sort { $best{$b} <=> $best{$a} || $a cmp $b } keys %best), "\n";
    }
    print "LD_SOURCE\tPLINK2_1KG_DIRECT\nLD_MIN_R2\t$threshold\n";
    print "LD_ESTIMABILITY\t", (keys(%best) ? 'ESTIMABLE' : 'NOT_ESTIMABLE'), "\n";
    if ($out) {
        my $population_label = population_label($populations, $keep);
        print {$out} join("\t", $ref, $ref, $population_label, 1, 'PLINK2_1KG_DIRECT', '1000_GENOMES_PHASE_3', $reference_build, ($phased ? 'PLINK2_R2_PHASED' : 'PLINK2_R2_UNPHASED')), "\n"
            if keys %best;
        for my $proxy (keys %best) {
            print {$out} join("\t", $ref, $proxy, $population_label, $best{$proxy}, 'PLINK2_1KG_DIRECT', '1000_GENOMES_PHASE_3', $reference_build, ($phased ? 'PLINK2_R2_PHASED' : 'PLINK2_R2_UNPHASED')), "\n";
        }
        close $out;
    }
}

sub population_label {
    my ($requested, $keep_file) = @_;
    if (defined($requested) && length($requested)) {
        my @pops = grep { length } map {
            my $x = uc($_);
            $x =~ s/^\s+|\s+$//g;
            $x;
        } split /,/, $requested;
        return join('+', @pops) if @pops;
    }
    return 'CUSTOM' if defined($keep_file) && length($keep_file);
    return 'ALL';
}

sub value {
    my ($f, $idx, @names) = @_;
    for my $name (@names) { return $f->[$idx->{$name}] if exists $idx->{$name}; }
    return undef;
}

sub native_path {
    my ($path) = @_;
    return $path unless defined($path) && length($path);
    my $abs = abs_path($path) || File::Spec->rel2abs($path);
    if ($abs =~ m{^/mnt/([A-Za-z])/(.*)$}) {
        my ($drive, $rest) = (uc($1), $2);
        $rest =~ s{/}{\\}g;
        return "$drive:\\$rest";
    }
    if ($^O =~ /cygwin/i && $abs =~ m{^/}) {
        if (open my $cygpath, '-|', 'cygpath', '-w', $abs) {
            my $win = <$cygpath> // '';
            close $cygpath;
            $win =~ s/[\r\n]+\z//;
            return $win if length $win;
        }
    }
    return $abs;
}

sub usage {
    return <<'USAGE';
Usage: resolve_plink2_local_ld.pl --query-snp rs123 (--pfile PREFIX | --bfile PREFIX) [options]
  --plink2 EXE       PLINK 2 executable (default: PLINK2 or plink2)
  --chr CHR --from-bp N --to-bp N   Optional local genomic interval
  --window-kb N      Maximum LD distance (default 1000)
  --min-r2 N         Minimum reported r2 (default 0; retain every in-window pair)
  --keep FILE        Optional PLINK sample keep file
  --populations LIST Restrict samples to 1000G superpopulations (EUR,AMR,AFR,EAS)
  --unphased         Use dosage-correlation r2 instead of phased haplotype r2
  --quiet            Suppress the potentially long LD_SNPS/LD_R2_PAIRS lines
  --output FILE      Write normalized LD cache rows
USAGE
}
