#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Cwd qw(abs_path);
use File::Path qw(make_path);
use File::Spec;
use Getopt::Long qw(GetOptions);
use Digest::MD5;
use JSON::PP;

# Public PGC schizophrenia release, European female and male autosomes + X.
# The expected MD5s are published by https://api.figshare.com/v2/articles/19426775.
my @files = (
 ['daner_PGC_SCZ_w3_75_0618a_eur_female.gz',41008148,'962c2bd94585b6af5cd9ecce1aa6d834'],
 ['daner_scz_w3_HRC_chrX_eur_fem_deduped_0518e.gz',52744046,'42d551bc0bf06f9323ae19ce69b119ce'],
 ['daner_PGC_SCZ_w3_75_0618a_eur_male.gz',41008139,'564694661075c1a805aa10586f2fc1f6'],
 ['daner_scz_w3_HRC_chrX_eur_mal_deduped_0518e.gz',52744043,'0a3fed754eaf7b5857d7214fb708aaf5'],
);
my ($output, $input, $gtf_cache, $plink2_1kg_pfile, $plink2, $help);
GetOptions('output-dir=s'=>\$output, 'input-dir=s'=>\$input,
           'gtf-cache-dir=s'=>\$gtf_cache,
           'plink2-1kg-pfile=s'=>\$plink2_1kg_pfile,
           'plink2=s'=>\$plink2, 'help'=>\$help) or die "Invalid options\n";
if ($help || !$output) {
 print "Usage: perl install/prepare_public_sex_gwas.pl --output-dir DIR [--input-dir EXISTING_GWAS_DIR] [--gtf-cache-dir DIR] [--plink2-1kg-pfile PREFIX --plink2 EXE]\n";
 exit($help ? 0 : 2);
}
$|=1;
my $reuse = defined $input;
make_path($output);
$output = abs_path($output);
$input ||= "$output/input";
make_path($input);
$input = abs_path($input);
my @provenance;
for my $f (@files) {
 my ($name,$id,$expected)=@$f;
 my $path="$input/$name";
 my $url="https://ndownloader.figshare.com/files/$id";
 my $downloaded=0;
 if (!-f $path) {
 die "Missing source file: $path\n" if $reuse;
 $path.='.part';
 print "Downloading $url\n";
  my @curl = ('curl','--fail','--location','--retry','3');
  push @curl, '--insecure'
    if ($ENV{PIPELINE_CURL_INSECURE} // '') =~ /^(?:1|true|yes|y|on)$/i;
  system(@curl,'--output',$path,$url)==0
    or die "Download failed: $url\n";
  $downloaded=1;
 }
 open my $fh,'<:raw',$path or die "Read $path: $!\n";
 my $md5=Digest::MD5->new->addfile($fh)->hexdigest;
 close $fh or die "Close $path: $!\n";
 die "Checksum mismatch: $path (expected $expected, got $md5)\n" unless $md5 eq $expected;
 rename($path,"$input/$name") or die "Rename download: $!\n" if $downloaded;
 push @provenance,{name=>$name,url=>$url,md5=>$md5,bytes=>-s "$input/$name"};
 print "Checksum PASS: $name\n";
}
my $spec={
 source_mode=>'raw_pgc_vcf_sumstats',project_tag=>'PUBLIC_SCZ_EUR_SEX',
 artifact_stem=>'public_scz_eur_sex',input_dir=>$input,output_dir=>$output,
 workdir=>abs_path("$Bin/.."),configs_dir=>"$output/configs",cygwin_bash=>'/bin/bash',reference_build=>'hg19',
 exclude_strand_ambiguous=>1,max_eaf_abs_diff=>0.2,threshold=>0.05,rho=>0,
 manhattan_differential_p_mode=>'raw',top_hit_focus_prefix=>'EUR',
 top_hit_maf_threshold=>0.01,top_hit_selection_method=>'ld',
 top_hit_ld_source=>'PLINK2_1KG',top_hit_ld_populations=>'EUR',
 local_ld_population=>'EUR',local_ld_display_mode=>'heatmap',local_ld_r2_threshold=>0.1,
 top_hit_ld_query_failure_action=>'KEEP',top_hit_max_loci=>3,
 local_window_bp=>500000,local_gtf_window_bp=>500000,
 open_result=>0,clean_oda_input=>1,keep_remote_plot_data=>0,
 groups=>[{tag=>'EUR_FEMALE',files=>[$files[0][0],$files[1][0]]},
           {tag=>'EUR_MALE',files=>[$files[2][0],$files[3][0]]}],
 pairs=>[{pair_tag=>'EUR_FEMALE_vs_MALE',group1=>'EUR_FEMALE',
          group2=>'EUR_MALE',prefix=>'EUR',label=>'EUR'}],
};
if (defined($plink2_1kg_pfile) || defined($plink2)) {
 die "Provide both --plink2-1kg-pfile and --plink2\n"
   unless defined($plink2_1kg_pfile) && defined($plink2);
 my $pgen_abs=abs_path("$plink2_1kg_pfile.pgen");
 my $prefix=defined($pgen_abs) ? $pgen_abs : $plink2_1kg_pfile;
 $prefix =~ s/\.pgen\z//;
 die "Missing 1000 Genomes Phase 3 PLINK2 files for prefix $prefix\n"
   unless -s "$prefix.pgen" && (-s "$prefix.pvar" || -s "$prefix.pvar.zst") && -s "$prefix.psam";
 my $plink_exe=abs_path($plink2) || $plink2;
 die "PLINK2 executable not found: $plink2\n" unless -f $plink_exe;
 $spec->{top_hit_ld_pfile}=$prefix;
 $spec->{top_hit_ld_plink2}=$plink_exe;
}
if (defined $gtf_cache) {
 die "GTF cache does not exist: $gtf_cache\n" unless -d $gtf_cache;
 $spec->{gtf_cache_dir}=abs_path($gtf_cache);
}
write_json("$output/spec.json",$spec);
write_json("$output/source_manifest.json",{
 study=>'Trubetskoy et al. Nature 2022, doi:10.1038/s41586-022-04434-5',
 source=>'https://figshare.com/articles/dataset/scz2022/19426775',
 build=>'GRCh37/hg19',comparison=>'European female minus European male',
 covariance_assumption=>'rho=0; disjoint sex strata',files=>\@provenance,
});
print "Spec: $output/spec.json\n";
sub write_json {
 my ($path,$data)=@_;
 open my $fh,'>',$path or die "Write $path: $!\n";
 print {$fh} JSON::PP->new->canonical->pretty->encode($data);
 close $fh or die "Close $path: $!\n";
}
