#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);
use JSON::PP;
use POSIX qw(strftime);

my ($out,$input,$cache,$plink2_1kg_pfile,$plink2,$help);
my $phase='all';
my $backend='gnuplot';
GetOptions('output-dir=s'=>\$out,'input-dir=s'=>\$input,'gtf-cache-dir=s'=>\$cache,
 'plink2-1kg-pfile=s'=>\$plink2_1kg_pfile,'plink2=s'=>\$plink2,
 'phase=s'=>\$phase,'backend=s'=>\$backend,'help'=>\$help) or die "Invalid options\n";
if ($help || !$out) {
 print <<'USAGE';
Usage: perl install/run_public_sex_gwas.pl --output-dir DIR [options]
  --input-dir DIR       Reuse four public source files after checksum checks
  --gtf-cache-dir DIR   Reuse an existing GENCODE hg19 cache
  --plink2-1kg-pfile P  1000 Genomes Phase 3 GRCh37 PLINK2 prefix
  --plink2 EXE          PLINK2 executable used for phased LD r2
  --phase NAME         all, prepare, preprocess, validate, plots, images
  --backend NAME       gnuplot (default), sas, both
Activate install/common.sh before running. SAS requires configured ODA access.
All phases use real female/male PGC schizophrenia GWAS, autosomes and chrX.
The full run requires several GB of space and can take hours on a slow host.
USAGE
 exit($help ? 0 : 2);
}
die "Invalid phase\n" unless $phase=~/^(?:all|prepare|preprocess|validate|plots|images)$/;
die "Invalid backend\n" unless $backend=~/^(?:gnuplot|sas|both)$/;
make_path($out); $out=abs_path($out);
$input=abs_path($input) if defined $input;
$cache=abs_path($cache) if defined $cache;
chdir "$Bin/.." or die $!;
$|=1;
my @steps;
my $spec="$out/spec.json";
if ($phase eq 'all' || $phase eq 'prepare') {
 my @cmd=($^X,"$Bin/prepare_public_sex_gwas.pl",'--output-dir',$out);
 push @cmd,('--input-dir',$input) if defined $input;
 push @cmd,('--gtf-cache-dir',$cache) if defined $cache;
 push @cmd,('--plink2-1kg-pfile',$plink2_1kg_pfile,'--plink2',$plink2)
   if defined($plink2_1kg_pfile) || defined($plink2);
 run('prepare',@cmd);
}
if ($phase eq 'all' || $phase eq 'preprocess') {
 run('preprocess',$^X,'auto_prepare_and_run_diff_gwas.pl','--spec',$spec,'--skip-plots');
}
if ($phase eq 'all' || $phase eq 'validate') {
 run('numeric_validation',$^X,"$Bin/validate_public_sex_gwas.pl",'--output-dir',$out);
}
if ($phase eq 'all' || $phase eq 'plots') {
 open my $fh,'<',"$out/targets.txt" or die "Run --phase validate first: $!\n";
 my $targets=<$fh>; close $fh; chomp $targets;
 die "Invalid target SNP list\n" unless $targets=~/^rs\d+(?:,rs\d+)*$/;
 if ($backend eq 'gnuplot' || $backend eq 'both') {
  run('gnuplot_inquiry',$^X,'auto_prepare_and_run_diff_gwas_with_gunplot.pl',
   '--spec',$spec,'--plots','manhattan,local_manhattan,local_gtf,forest',
   '--target-snps',$targets,'--no-remove-X-chr');
  verify_images('GUNPLOT');
 }
 if ($backend eq 'sas' || $backend eq 'both') {
  run('sas_login',$^X,'run_sas_codes_or_script_in_ODA.pl','--check-sas-oda-login-only');
  run('sas_inquiry',$^X,'auto_prepare_and_run_diff_gwas.pl',
   '--spec',$spec,'--plots','manhattan,local_manhattan,local_gtf,forest',
   '--target-snps',$targets,'--from-step','plot_manhattan','--force',
   '--no-gnuplot-fallback-on-sas-space','--no-gnuplot-fallback-on-sas-failure');
  verify_images('SAS');
 }
}
sub verify_images {
 my ($backend_name)=@_;
 require GD;
 open my $sf,'<',$spec or die $!;
 my $cfg=decode_json(do {local $/;<$sf>}); close $sf;
 my $prefix=$cfg->{project_tag}.'_'.$backend_name.'_';
 my (@images,%families);
 for my $folder ($out,abs_path('.')) {
  opendir my $dh,$folder or die $!;
  my @names=grep {index($_,$prefix)==0 && /\.png$/i} readdir $dh;
  closedir $dh;
  for my $name (@names) {
   my $path="$folder/$name";
   open my $fh,'<:raw',$path or die $!;
   my $im=GD::Image->new($fh) or die "Invalid PNG: $path\n";
   close $fh;
   die "Unexpectedly small PNG: $path\n" if $im->width<100 || $im->height<100;
   my $family=$name=~/forest/i?'forest':$name=~/gtf/i?'local_gtf':$name=~/local/i?'local_manhattan':'manhattan';
   $families{$family}++;
   push @images,{path=>$path,width=>$im->width,height=>$im->height,bytes=>-s $path,family=>$family};
  }
 }
 for (qw(manhattan local_manhattan local_gtf forest)) {
  die "Missing $backend_name $_ PNG output\n" unless $families{$_};
 }
 open my $vf,'>',"$out/image_validation_".lc($backend_name).'.json' or die $!;
 print {$vf} JSON::PP->new->canonical->pretty->encode({status=>'PASS',images=>\@images});
 close $vf or die $!;
 print "PASS: decoded $backend_name PNGs for all four plot families\n";
}
if ($phase eq 'images') {
 verify_images('GUNPLOT') if $backend eq 'gnuplot' || $backend eq 'both';
 verify_images('SAS') if $backend eq 'sas' || $backend eq 'both';
}
if ($phase eq 'all' || $phase eq 'plots' || $phase eq 'images') {
 system($^X,"$Bin/build_public_gwas_gallery.pl",'--output-dir',$out)==0
   or die "Building results gallery failed\n";
}
print "Requested phase completed. Reports: $out\n";
sub run {
 my ($name,@cmd)=@_;
 print "\n[$name] ",join(' ',@cmd),"\n";
 my $start=time;
 my $rc=system(@cmd);
 push @steps,{name=>$name,command=>\@cmd,seconds=>time-$start,raw_status=>$rc,
  status=>$rc==0?'PASS':'FAIL',finished_utc=>strftime('%Y-%m-%dT%H:%M:%SZ',gmtime)};
 open my $fh,'>',"$out/test_run_$phase.json" or die $!;
 print {$fh} JSON::PP->new->canonical->pretty->encode({steps=>\@steps});
 close $fh or die $!;
 die "$name failed (status $rc); see the command output above\n" if $rc;
}
