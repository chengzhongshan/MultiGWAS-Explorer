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
if (!defined($plink2_1kg_pfile) && !defined($plink2)) {
 my $candidate_pfile=abs_path('cache/plink2_1kg_phase3/all_phase3')
   || 'cache/plink2_1kg_phase3/all_phase3';
 my @plink_names=$^O eq 'MSWin32' ? qw(plink2.exe plink2) : qw(plink2 plink2.exe);
 my ($candidate_plink)=map { abs_path("cache/plink2_bin/$_") }
   grep { -f "cache/plink2_bin/$_" } @plink_names;
 if (-s "$candidate_pfile.pgen"
     && (-s "$candidate_pfile.pvar" || -s "$candidate_pfile.pvar.zst")
     && -s "$candidate_pfile.psam" && defined($candidate_plink)) {
  ($plink2_1kg_pfile,$plink2)=($candidate_pfile,$candidate_plink);
  print "Using repository-local PLINK2 Phase 3 reference: $candidate_pfile\n";
 }
}
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
 validate_plink_reference($spec);
 if ($backend eq 'gnuplot' || $backend eq 'both') {
  run('gnuplot_inquiry',$^X,'auto_prepare_and_run_diff_gwas_with_gunplot.pl',
   '--spec',$spec,'--plots','manhattan,local_manhattan,local_gtf,forest',
   '--target-snps',$targets,'--no-remove-X-chr','--ld-display-mode','heatmap');
  verify_images('GUNPLOT');
  verify_gnuplot_signed_ld($targets);
 }
 if ($backend eq 'sas' || $backend eq 'both') {
  run('sas_login',$^X,'run_sas_codes_or_script_in_ODA.pl','--check-sas-oda-login-only');
  run('sas_inquiry',$^X,'auto_prepare_and_run_diff_gwas.pl',
   '--spec',$spec,'--plots','manhattan,local_manhattan,local_gtf,forest',
   '--target-snps',$targets,'--from-step','plot_manhattan','--force',
   '--local-ld-display-mode','heatmap','--local-ld-population','EUR',
   '--local-ld-r2-threshold','0.1',
   '--no-gnuplot-fallback-on-sas-space','--no-gnuplot-fallback-on-sas-failure');
  verify_images('SAS');
  verify_sas_target_gtf($targets);
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
sub validate_plink_reference {
 my ($spec_path)=@_;
 open my $fh,'<',$spec_path or die $!;
 my $cfg=decode_json(do {local $/;<$fh>}); close $fh;
 my $prefix=$cfg->{top_hit_ld_pfile}//'';
 my $exe=$cfg->{top_hit_ld_plink2}//'';
 die "The public GTF example requires a configured PLINK2 Phase 3 reference. "
   . "Rerun --phase prepare with --plink2-1kg-pfile and --plink2.\n"
  unless length($prefix) && -s "$prefix.pgen"
    && (-s "$prefix.pvar" || -s "$prefix.pvar.zst") && -s "$prefix.psam"
    && length($exe) && -f $exe;
}
sub verify_gnuplot_signed_ld {
 my ($targets)=@_;
 my @wanted=split /,/,$targets;
 opendir my $dh,$out or die "Cannot inspect $out: $!\n";
 my @files=grep {/GUNPLOT_local_top_hits_with_gtf_.*\.manifest\.tsv\z/} readdir $dh;
 closedir $dh;
 for my $snp (@wanted) {
  my ($file)=grep {/\Q$snp\E\.manifest\.tsv\z/} @files;
  die "Missing gnuplot local-GTF manifest for $snp\n" unless $file;
  open my $fh,'<',"$out/$file" or die $!;
  my %metric;
  while (<$fh>) { chomp; s/\r\z//; my ($k,$v)=split /\t/,$_,2; $metric{$k}=$v if defined $v; }
  close $fh;
  die "$snp did not use heatmap LD display\n" unless ($metric{ld_display_mode}//'') eq 'heatmap';
  die "$snp did not render r2 x sign(Z)\n" unless ($metric{signed_r2_coloring}//0)==1;
  die "$snp has no PLINK2 LD proxies in the plotted locus\n" unless ($metric{ld_r2_points}//0)>1;
  die "$snp LD values did not come from the PLINK2 Phase 3 cache\n"
   unless ($metric{ld_source_file}//'') =~ /\.plink2_1kg_phase3\.tsv\z/;
  my $source=$metric{ld_source_file};
  die "$snp PLINK2 Phase 3 cache is missing\n" unless -s $source;
  open my $lf,'<',$source or die $!;
  my @header=split /\t/,scalar(<$lf>),-1;
  my @first=split /\t/,scalar(<$lf>),-1;
  close $lf;
  s/[\r\n]+\z// for @header,@first;
  my %idx=map {$header[$_]=>$_} 0..$#header;
  die "$snp LD cache has an incomplete provenance header\n"
   unless !grep {!exists $idx{$_}} qw(ld_population source reference_build ld_method);
  die "$snp LD cache is not phased EUR Phase 3 GRCh37 output\n"
   unless ($first[$idx{ld_population}]//'') eq 'EUR'
    && ($first[$idx{source}]//'') eq 'PLINK2_1KG_DIRECT'
    && ($first[$idx{reference_build}]//'') eq 'GRCh37_hg19'
    && ($first[$idx{ld_method}]//'') eq 'PLINK2_R2_PHASED';
 }
 print "PASS: gnuplot local-GTF panels use PLINK2 Phase 3 r2 x sign(Z)\n";
}
sub verify_sas_target_gtf {
 my ($targets)=@_;
 require GD;
 open my $sf,'<',$spec or die $!;
 my $cfg=decode_json(do {local $/;<$sf>}); close $sf;
 my $prefix=$cfg->{project_tag}.'_SAS_local_top_hits_with_gtf_';
 my @folders=($out,abs_path('.'));
 for my $snp (split /,/,$targets) {
  (my $safe_snp=$snp)=~s/[^A-Za-z0-9._-]/_/g;
  my $name=$prefix.$safe_snp.'.png';
  my ($path)=grep {-s $_} map {"$_/$name"} @folders;
  die "Missing target-specific SAS local-GTF PNG for $snp: $name\n" unless $path;
  open my $fh,'<:raw',$path or die $!;
  my $im=GD::Image->new($fh) or die "Invalid target-specific SAS PNG: $path\n";
  close $fh;
  die "Unexpectedly small target-specific SAS PNG: $path\n"
   if $im->width<100 || $im->height<100;
 }
 print "PASS: SAS produced one target-specific local-GTF panel per SNP\n";
}
if ($phase eq 'images') {
 if ($backend eq 'gnuplot' || $backend eq 'both') {
  verify_images('GUNPLOT');
  open my $tf,'<',"$out/targets.txt" or die "Run --phase validate first: $!\n";
  my $targets=<$tf>; close $tf; chomp $targets;
  verify_gnuplot_signed_ld($targets);
 }
 if ($backend eq 'sas' || $backend eq 'both') {
  verify_images('SAS');
  open my $tf,'<',"$out/targets.txt" or die "Run --phase validate first: $!\n";
  my $targets=<$tf>; close $tf; chomp $targets;
  verify_sas_target_gtf($targets);
 }
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
