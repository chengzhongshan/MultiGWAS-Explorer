#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw($GunzipError);
use File::Path qw(make_path);

my ($input,$outdir,$prefix,$nmax)=('','benchmark/fixture','ALL',100_000);
GetOptions('input=s'=>\$input,'output-dir=s'=>\$outdir,'prefix=s'=>\$prefix,'rows=i'=>\$nmax) or die usage();
die usage() unless $input && $nmax>0; make_path($outdir);
my $fh=$input=~/\.gz$/i ? IO::Uncompress::Gunzip->new($input, MultiStream=>1)||die $GunzipError : do{open(my$x,'<',$input)or die$!;$x};
my $h=<$fh>//die"Empty input\n"; chomp$h; $h=~s/^#//; my@c=split/\t/,$h,-1; my%i=map{$c[$_]=>$_}0..$#c;
for(qw(CHR BP A1 A2 SNP GWAS_TAG BETA SE FRQ_A FRQ_U NEFF)){die"Missing $_\n"unless exists$i{$_}}
my($g1,$g2)=("${prefix}_FEMALE","${prefix}_MALE");
open my$f,'>',"$outdir/female.gwama.tsv"or die$!; open my$m,'>',"$outdir/male.gwama.tsv"or die$!;
open my$e,'>',"$outdir/easystrata.tsv"or die$!; open my$x,'>',"$outdir/expected.tsv"or die$!;
open my$w,'>',"$outdir/multigwas.tsv"or die$!;
print $_ "MARKERNAME\tSTRAND\tCHR\tPOS\tIMPUTED\tEA\tNEA\tEAF\tN\tBETA\tSE\n" for($f,$m);
print $e "SNP\tCHR\tBP\tA1\tA2\tBETA_FEMALE\tSE_FEMALE\tBETA_MALE\tSE_MALE\n";
print $x "SNP\tCHR\tBP\tDIFF_BETA\tDIFF_SE\tDIFF_Z\tDIFF_P\tQ\n";
print $w "CHR\tBP\tA1\tA2\tSNP\tALL_GROUP1_BETA\tALL_GROUP2_BETA\tALL_DIFF_BETA\tALL_GROUP1_SE\tALL_GROUP2_SE\tALL_DIFF_SE\tALL_DIFF_P\n";
my $key = '';
my %b;
my $n = 0;
while(my$line=<$fh>){chomp$line;my@v=split/\t/,$line,-1;my$k=join("\t",map{$v[$i{$_}]//''}qw(CHR BP A1 A2 SNP));
 if(length$key&&$k ne$key){emit();%b=();last if$n>=$nmax}$key=$k;my$t=$v[$i{GWAS_TAG}]//'';$b{$t}=[@v]if$t eq$g1||$t eq$g2;}
emit() if $n < $nmax && %b; close $_ for($fh,$f,$m,$e,$x,$w);
open my$l,'>',"$outdir/gwama.in"or die$!; print {$l} "female.gwama.tsv F\nmale.gwama.tsv M\n"; close $l;
open my$c,'>',"$outdir/easystrata.ecf"or die$!; print {$c} <<'ECF';
# EasyStrata 8.6 has an upstream empty-comment-index parser bug; keep this row.
DEFINE --pathOut ./
       --strSeparator TAB
EASYIN --fileIn easystrata.tsv --fileInShortName easystrata
START EASYX
CALCPDIFF --acolBETAs BETA_FEMALE;BETA_MALE
          --acolSEs SE_FEMALE;SE_MALE
          --colOutPdiff EASYSTRATA_PDIFF
WRITE --strMode txt
      --strPrefix EasyStrata.
      --strSuffix .results
      --strSep TAB
STOP EASYX
ECF
close$c;print"Exported $n matched $prefix pairs to $outdir\n";
sub emit{return unless$b{$g1}&&$b{$g2};my($a,$d)=@b{$g1,$g2};my($b1,$b2,$s1,$s2)=map{num($_)}($a->[$i{BETA}],$d->[$i{BETA}],$a->[$i{SE}],$d->[$i{SE}]);return unless defined$b1&&defined$b2&&defined$s1&&defined$s2&&$s1>0&&$s2>0;
 my($chr,$bp,$a1,$a2,$snp)=map{$a->[$i{$_}]//''}qw(CHR BP A1 A2 SNP);my$ef1=meanfreq($a->[$i{FRQ_A}],$a->[$i{FRQ_U}]);my$ef2=meanfreq($d->[$i{FRQ_A}],$d->[$i{FRQ_U}]);my$n1=num($a->[$i{NEFF}])//1;my$n2=num($d->[$i{NEFF}])//1;
 print$f join("\t",$snp,'+',$chr,$bp,1,$a1,$a2,$ef1,$n1,$b1,$s1),"\n";print$m join("\t",$snp,'+',$chr,$bp,1,$a1,$a2,$ef2,$n2,$b2,$s2),"\n";print$e join("\t",$snp,$chr,$bp,$a1,$a2,$b1,$s1,$b2,$s2),"\n";my$db=$b1-$b2;my$ds=sqrt($s1*$s1+$s2*$s2);my$z=$db/$ds;my$p=erfc(abs($z)/sqrt(2));print$x join("\t",$snp,$chr,$bp,$db,$ds,$z,$p,$z*$z),"\n";print$w join("\t",$chr,$bp,$a1,$a2,$snp,$b1,$b2,$db,$s1,$s2,$ds,$p),"\n";$n++;}
sub num{defined$_[0]&&$_[0]=~/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i?0+$_[0]:undef}
sub meanfreq{my@q=grep{defined}map{num($_)}@_;return@q?eval(join('+',@q))/@q:0.5}
sub erfc{my($z)=@_;my$t=1/(1+.5*$z);$t*exp(-$z*$z-1.26551223+$t*(1.00002368+$t*(.37409196+$t*(.09678418+$t*(-.18628806+$t*(.27886807+$t*(-1.13520398+$t*(1.48851587+$t*(-.82215223+$t*.17087277)))))))))}
sub usage{"Usage: perl export_official_benchmark_fixture.pl --input MERGED_LONG.tsv.gz [--output-dir DIR] [--prefix ALL] [--rows 100000]\n"}
