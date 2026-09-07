#!/usr/bin/env perl
use strict;
use warnings;

use File::Temp qw(tempdir);
use Getopt::Long qw(GetOptions);
use Text::ParseWords qw(shellwords);

my ($query, $pfile, $plink2, $chr, $from_bp, $to_bp, $keep, $output);
my $min_r2 = 0.8;
my $window_kb = 1000;
my $phased = 1;

GetOptions(
    'query-snp=s' => \$query,
    'pfile=s' => \$pfile,
    'plink2=s' => \$plink2,
    'chr=s' => \$chr,
    'from-bp=i' => \$from_bp,
    'to-bp=i' => \$to_bp,
    'window-kb=f' => \$window_kb,
    'min-r2=f' => \$min_r2,
    'keep=s' => \$keep,
    'unphased' => sub { $phased = 0 },
    'output=s' => \$output,
) or die usage();

die usage() unless defined($query) && defined($pfile);
die "--min-r2 must be between 0 and 1\n" unless $min_r2 > 0 && $min_r2 <= 1;
die "--window-kb must be positive\n" unless $window_kb > 0;
$plink2 ||= $ENV{PLINK2} || 'plink2';

my $tmp = tempdir('plink2_local_ld.XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $prefix = "$tmp/ld";
my @cmd = (shellwords($plink2));
push @cmd, '--pfile', $pfile, 'vzs' if -e "$pfile.pvar.zst";
push @cmd, '--pfile', $pfile unless -e "$pfile.pvar.zst";
push @cmd, '--chr', $chr if defined($chr) && length($chr);
push @cmd, '--from-bp', $from_bp if defined $from_bp;
push @cmd, '--to-bp', $to_bp if defined $to_bp;
push @cmd, '--keep', $keep if defined($keep) && length($keep);
push @cmd, '--ld-snp', $query;
push @cmd, ($phased ? '--r2-phased' : '--r2-unphased');
push @cmd, '--ld-window-kb', $window_kb;
push @cmd, '--ld-window-r2', $min_r2;
push @cmd, '--allow-ambiguous-allele', '--out', $prefix;

system(@cmd);
die "PLINK2 LD calculation failed (exit $?): @cmd\n" if $? != 0;
my ($report) = grep { -f $_ } ("$prefix.vcor", "$prefix.vcor.zst");
die "PLINK2 did not produce a vcor report\n" unless $report;
if ($report =~ /\.zst\z/) {
    my $zcat = $ENV{ZSTD} || 'zstd';
    open my $in, '-|', $zcat, '-dc', $report or die "Cannot decompress $report: $!\n";
    parse_report($in, $output, $query, $min_r2);
    close $in;
} else {
    open my $in, '<', $report or die "Cannot read $report: $!\n";
    parse_report($in, $output, $query, $min_r2);
    close $in;
}

sub parse_report {
    my ($fh, $out_path, $ref, $threshold) = @_;
    my ($header, %idx, %best);
    while (my $line = <$fh>) {
        chomp $line;
        next unless length $line;
        if (!$header) {
            $header = 1;
            my @h = split /\t/, $line, -1;
            %idx = map { lc($h[$_]) =~ s/^#//r => $_ } 0 .. $#h;
            next;
        }
        my @f = split /\t/, $line, -1;
        my $a = value(\@f, \%idx, qw(id_a variant_id_a));
        my $b = value(\@f, \%idx, qw(id_b variant_id_b));
        my $r2 = value(\@f, \%idx, qw(r2));
        next unless defined($a) && defined($b) && defined($r2);
        next unless $a =~ /^rs\d+$/i && $b =~ /^rs\d+$/i;
        next unless $r2 =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i;
        next unless $r2 >= $threshold;
        my $proxy = lc($a) eq lc($ref) ? $b : (lc($b) eq lc($ref) ? $a : undef);
        next unless defined $proxy;
        $best{$proxy} = 0 + $r2 if !exists($best{$proxy}) || $r2 > $best{$proxy};
    }
    open my $out, '>', $out_path or die "Cannot write $out_path: $!\n" if defined $out_path;
    print {$out} join("\t", qw(query_snp proxy_snp ld_population proxy_r2 source)), "\n" if $out;
    print "LD_SNPS\t", join(',', sort { $best{$b} <=> $best{$a} || $a cmp $b } keys %best), "\n";
    print "LD_R2_PAIRS\t", join(',', map { $_ . ':' . sprintf('%.6g', $best{$_}) } sort { $best{$b} <=> $best{$a} || $a cmp $b } keys %best), "\n";
    print "LD_SOURCE\tPLINK2_LOCAL\nLD_MIN_R2\t$threshold\n";
    if ($out) {
        for my $proxy (keys %best) {
            print {$out} join("\t", $ref, $proxy, '1KG', $best{$proxy}, 'PLINK2_LOCAL'), "\n";
        }
        close $out;
    }
}

sub value {
    my ($f, $idx, @names) = @_;
    for my $name (@names) { return $f->[$idx->{$name}] if exists $idx->{$name}; }
    return undef;
}

sub usage {
    return <<'USAGE';
Usage: resolve_plink2_local_ld.pl --query-snp rs123 --pfile /path/chr2_phase3 [options]
  --plink2 EXE       PLINK 2 executable (default: PLINK2 or plink2)
  --chr CHR --from-bp N --to-bp N   Optional local genomic interval
  --window-kb N      Maximum LD distance (default 1000)
  --min-r2 N         Minimum reported r2 (default 0.8)
  --keep FILE        Optional PLINK sample keep file
  --unphased         Use dosage-correlation r2 instead of phased haplotype r2
  --output FILE      Write normalized LD cache rows
USAGE
}
