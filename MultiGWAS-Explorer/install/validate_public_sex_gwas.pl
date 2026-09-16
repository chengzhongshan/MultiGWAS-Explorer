#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw($GunzipError);
use JSON::PP;
use POSIX qw(erfc isfinite);

my $dir;
GetOptions('output-dir=s'=>\$dir) && $dir
 or die "Usage: perl install/validate_public_sex_gwas.pl --output-dir DIR\n";
$|=1;
my $path="$dir/public_scz_eur_sex.tsv.gz";
my $fh=IO::Uncompress::Gunzip->new($path, MultiStream=>1)
 or die "Read $path: $GunzipError\n";
my $header=<$fh> // die "Empty differential file\n";
chomp $header;
$header =~ s/\r$//;
my @cols=split /\t/,$header;
my %idx=map {$cols[$_]=>$_} 0..$#cols;
for (qw(CHR BP SNP GROUP1_BETA GROUP2_BETA GROUP1_SE GROUP2_SE DIFF_BETA DIFF_SE DIFF_Z DIFF_P GROUP1_P GROUP2_P STRAND_AMBIGUOUS EAF_ABS_DIFF)) {
 die "Missing column $_\n" unless exists $idx{$_};
}
my ($n,%chr,%best);
$n=0;
while (my $line=<$fh>) {
 chomp $line;
 my @f=split /\t/,$line,-1;
 die "Wrong field count at row $n\n" unless @f==@cols;
 my %r; @r{@cols}=@f;
 ++$n;
 for my $k (qw(GROUP1_BETA GROUP2_BETA GROUP1_SE GROUP2_SE DIFF_BETA DIFF_SE DIFF_Z DIFF_P)) {
  die "Invalid $k at $r{SNP}\n" unless $r{$k}=~/^[+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?$/i && isfinite(0+$r{$k});
 }
 die "Nonpositive SE at $r{SNP}\n" unless $r{GROUP1_SE}>0 && $r{GROUP2_SE}>0;
 my $beta=$r{GROUP1_BETA}-$r{GROUP2_BETA};
 my $se=sqrt($r{GROUP1_SE}**2+$r{GROUP2_SE}**2);
 my $z=$beta/$se;
 my $p=erfc(abs($z)/sqrt(2));
 near($r{DIFF_BETA},$beta,1e-8,$r{SNP},'beta');
 near($r{DIFF_SE},$se,1e-8,$r{SNP},'SE');
 near($r{DIFF_Z},$z,1e-7,$r{SNP},'Z');
 # The production code uses an erfc approximation and six-digit P output.
 die "P mismatch at $r{SNP}\n" if abs($r{DIFF_P}-$p)>2e-6*($p+1e-300);
 die "P outside [0,1] at $r{SNP}\n" unless $r{DIFF_P}>=0 && $r{DIFF_P}<=1;
 die "Ambiguous SNP survived: $r{SNP}\n" if $r{STRAND_AMBIGUOUS};
 die "EAF discordance survived: $r{SNP}\n" if length($r{EAF_ABS_DIFF}) && $r{EAF_ABS_DIFF}>0.2000001;
 $chr{$r{CHR}}++;
 if ($r{SNP}=~/^rs\d+$/) {
  my $key=($r{CHR}=~/^(?:23|X)$/i) ? 'chrX' : 'differential';
  $best{$key}={snp=>$r{SNP},chr=>$r{CHR},bp=>$r{BP},p=>$r{DIFF_P}}
    if !exists($best{$key}) || $r{DIFF_P}<$best{$key}{p};
  my $common=$r{GROUP1_P}<$r{GROUP2_P} ? $r{GROUP1_P}:$r{GROUP2_P};
  $best{common}={snp=>$r{SNP},chr=>$r{CHR},bp=>$r{BP},p=>$common}
    if !exists($best{common}) || $common<$best{common}{p};
 }
 print "Validated $n differential rows\n" unless $n%1000000;
}
die "Gzip read error: ".$fh->error."\n" if $fh->error;
$fh->close or die "Gzip close failed\n";
die "No differential rows\n" unless $n;
for (1..22) { die "Missing autosome $_\n" unless $chr{$_} }
die "Missing chromosome X\n" unless $chr{23} || $chr{X};
my $report={status=>'PASS',rows_checked=>$n,chromosomes=>\%chr,
 independent_numeric_check=>'beta_F-beta_M, sqrt(SE_F^2+SE_M^2), POSIX::erfc',
 inquiry_targets=>\%best};
open my $out,'>',"$dir/numeric_validation.json" or die $!;
print {$out} JSON::PP->new->canonical->pretty->encode($report);
close $out or die $!;
open my $targets,'>',"$dir/targets.txt" or die $!;
my %seen;
print {$targets} join(',',grep {!$seen{$_}++} map {$best{$_}{snp}} qw(common differential chrX)),"\n";
close $targets or die $!;
print "PASS: $n differential rows; autosomes 1-22 and X; independent effect/Z/P checks\n";
sub near {
 my ($got,$want,$tolerance,$snp,$label)=@_;
 die "$label mismatch at $snp: $got versus $want\n" if abs($got-$want)>$tolerance*(1+abs($want));
}
