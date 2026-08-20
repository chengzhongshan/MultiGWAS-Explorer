#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw($GunzipError);
use Time::HiRes qw(time);

my ($input, $output, $prefix, $limit, $rho) = ('', 'benchmark_diff_methods.tsv', 'ALL', 0, 0);
GetOptions('input=s'=>\$input, 'output=s'=>\$output, 'prefix=s'=>\$prefix,
           'limit=i'=>\$limit, 'rho=f'=>\$rho) or die usage();
die usage() unless $input;
die "rho must be in (-1,1)\n" unless $rho > -1 && $rho < 1;
my $fh = $input =~ /\.gz$/i
  ? IO::Uncompress::Gunzip->new($input) || die "Cannot open $input: $GunzipError\n"
  : do { open(my $x, '<', $input) or die $!; $x };
my $h = <$fh> // die "Empty input\n"; chomp $h; $h =~ s/\r$//;
my @c = split /\t/, $h, -1; my %i = map { $c[$_] => $_ } 0..$#c;
my @need = map { "${prefix}_$_" } qw(GROUP1_BETA GROUP2_BETA GROUP1_SE GROUP2_SE DIFF_BETA DIFF_SE DIFF_P);
die "Missing $_\n" for grep { !exists $i{$_} } @need;
my ($n,$max_b,$max_se,$max_z2,$max_p,$start)=(0,0,0,0,0,time);
while(my $line=<$fh>){ chomp $line; my @v=split /\t/,$line,-1;
  my ($b1,$b2,$s1,$s2)=map { num($v[$i{"${prefix}_$_"}]) } qw(GROUP1_BETA GROUP2_BETA GROUP1_SE GROUP2_SE);
  next unless defined $b1 && defined $b2 && defined $s1 && defined $s2 && $s1>0 && $s2>0;
  my $db=$b1-$b2; my $ds=sqrt($s1*$s1+$s2*$s2-2*$rho*$s1*$s2); my $z=$db/$ds;
  my $p=erfc(abs($z)/sqrt(2));
  # EasyStrata CALCPDIFF option 1 uses this beta-difference test. For two
  # independent strata, GWAMA Cochran Q is algebraically z_diff squared.
  my $q=$z*$z;
  $max_b=max($max_b,abs($db-num($v[$i{"${prefix}_DIFF_BETA"}])));
  $max_se=max($max_se,abs($ds-num($v[$i{"${prefix}_DIFF_SE"}])));
  my $stored_p=num($v[$i{"${prefix}_DIFF_P"}]); $max_p=max($max_p,abs($p-$stored_p)) if defined $stored_p;
  $max_z2=max($max_z2,abs($q-$z*$z)); $n++; last if $limit && $n >= $limit;
}
close $fh; my $sec=time-$start;
open my $o,'>',$output or die $!;
print {$o} "METRIC\tVALUE\nrows_compared\t$n\nelapsed_seconds\t$sec\nrows_per_second\t",($sec?$n/$sec:0),"\n";
print {$o} "max_abs_diff_beta\t$max_b\nmax_abs_diff_se\t$max_se\nmax_abs_diff_p\t$max_p\nmax_abs_diff_gwama_q_vs_z2\t$max_z2\n";
print {$o} "easystrata_target\tCALCPDIFF option 1, covariance correction rho=$rho\n";
print {$o} "gwama_target\ttwo-stratum Cochran Q equivalence when rho=0\n";
close $o; print "Compared $n rows in ",sprintf('%.3f',$sec)," seconds; wrote $output\n";
sub num { defined $_[0] && $_[0]=~/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i ? 0+$_[0] : undef }
sub max { $_[0] > $_[1] ? $_[0] : $_[1] }
sub erfc { my($x)=@_; my$t=1/(1+.5*$x); return $t*exp(-$x*$x-1.26551223+$t*(1.00002368+$t*(.37409196+$t*(.09678418+$t*(-.18628806+$t*(.27886807+$t*(-1.13520398+$t*(1.48851587+$t*(-.82215223+$t*.17087277))))))))) }
sub usage { "Usage: perl benchmark_diff_methods.pl --input WIDE.tsv[.gz] [--prefix ALL] [--rho 0] [--limit N] [--output FILE]\n" }
