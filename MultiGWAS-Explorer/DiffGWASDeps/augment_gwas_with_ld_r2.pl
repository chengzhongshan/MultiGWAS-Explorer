#!/usr/bin/env perl
use strict; use warnings; use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError);
my ($input,$cache,$output,$reference);
GetOptions('input=s'=>\$input,'ld-cache=s'=>\$cache,'output=s'=>\$output,'reference-snp=s'=>\$reference) or die usage();
die usage() unless $input && $cache && $output;
open my $cfh,'<',$cache or die "Cannot read LD cache $cache: $!\n";
my $ch=<$cfh>; die "Empty LD cache\n" unless defined $ch; chomp $ch; $ch =~ s/\r$//; my @hc=split /\t/,$ch,-1; my %ci=map {lc($hc[$_])=>$_} 0..$#hc;
die "LD cache requires proxy_snp and proxy_r2 columns\n" unless exists $ci{proxy_snp} && exists $ci{proxy_r2};
my %r2; while(my $line=<$cfh>){chomp $line; $line=~s/\r$//; my @f=split /\t/,$line,-1;
    if (defined($reference) && length($reference) && exists($ci{query_snp})) {
        my $query=$f[$ci{query_snp}]//'';
        next unless lc($query) eq lc($reference);
    }
    my $s=$f[$ci{proxy_snp}]//' '; my $r=$f[$ci{proxy_r2}]//' '; next unless $s ne ' ' && $r =~ /^\d/; $r2{lc $s}=0+$r if !exists($r2{lc $s}) || $r>$r2{lc $s};} close $cfh;
$r2{lc $reference}=1 if defined $reference && length $reference;
my $in=$input=~/\.gz$/ ? IO::Uncompress::Gunzip->new($input) : do {open my $fh,'<',$input or die $!; $fh}; die "Cannot read input\n" unless $in;
my $out=$output=~/\.gz$/ ? IO::Compress::Gzip->new($output) : do {open my $fh,'>',$output or die $!; $fh}; die "Cannot write output\n" unless $out;
my $hdr=<$in>; die "Empty input\n" unless defined $hdr; chomp $hdr; $hdr=~s/\r$//; my @h=split /\t/,$hdr,-1; my %i=map {uc($h[$_])=>$_} 0..$#h; die "Input requires SNP column\n" unless exists $i{SNP}; my $has=exists $i{LD_R2}; push @h,'LD_R2' unless $has; print $out join("\t",@h),"\n";
while(my $line=<$in>){chomp $line; $line=~s/\r$//; my @f=split /\t/,$line,-1; my $v=$r2{lc($f[$i{SNP}]//'')}; if($has){$f[$i{LD_R2}]=defined($v)?sprintf('%.12g',$v):'';} else {push @f,defined($v)?sprintf('%.12g',$v):'';} print $out join("\t",@f),"\n";} close $in; close $out; print "AUGMENTED_ROWS\t",scalar(keys %r2),"\n";
sub usage { 'Usage: augment_gwas_with_ld_r2.pl --input DATA --ld-cache CACHE --output OUT [--reference-snp SNP]' }
