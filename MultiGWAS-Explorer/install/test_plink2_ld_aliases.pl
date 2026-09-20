#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use Test::More;

my $dir = tempdir(CLEANUP => 1);
my $fake = "$dir/fake_plink2.pl";
open my $fh, '>', $fake or die $!;
print {$fh} <<'FAKE';
#!/usr/bin/env perl
use strict;
use warnings;
my ($out);
for (my $i = 0; $i < @ARGV; $i++) {
    $out = $ARGV[$i + 1] if $ARGV[$i] eq '--out';
}
die "missing --out\n" unless defined $out;
if ($^O =~ /cygwin/i && $out =~ /^[A-Za-z]:[\\\/]/) {
    open my $cp, '-|', 'cygpath', '-u', $out or die $!;
    $out = <$cp> // '';
    close $cp;
    $out =~ s/[\r\n]+\z//;
}
open my $v, '>', "$out.vcor" or die $!;
print {$v} "#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tR2\n";
print {$v} "1\t100\trsQ;rsAlias\t1\t110\trsA;rsB\t0.8\n";
print {$v} "1\t100\trsQ\t1\t120\t.\t0.7\n";
print {$v} "1\t130\trsOther\t1\t140\trsUnused\t0.9\n";
close $v;
FAKE
close $fh;
chmod 0755, $fake or die $!;

my $output = "$dir/ld.tsv";
my @cmd = ($^X, "$Bin/../DiffGWASDeps/resolve_plink2_local_ld.pl",
    '--query-snp', 'rsQ', '--pfile', "$dir/ref", '--plink2', $fake,
    '--min-r2', '0.1', '--output', $output, '--quiet');
is(system(@cmd), 0, 'PLINK2 LD normalizer accepts the fixture report');
open $fh, '<', $output or die $!;
my @lines = <$fh>;
close $fh;
my $text = join('', @lines);
like($text, qr/^rsQ\trsA\t/m, 'first semicolon-delimited alias is emitted');
like($text, qr/^rsQ\trsB\t/m, 'second semicolon-delimited alias is emitted');
unlike($text, qr/\t\.\t/, 'missing-ID sentinel is discarded');
unlike($text, qr/rsA;rsB/, 'compound alias is not emitted as one unusable ID');
done_testing;
