#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use JSON::PP;
use Test::More;
my $dir=tempdir('multigwas config XXXXX',TMPDIR=>1,CLEANUP=>1);
my $source="$dir/existing.tsv";
open my $fh,'>',$source or die $!;
print {$fh} "CHR\tBP\tSNP\tGROUP1_BETA\tGROUP2_BETA\tDIFF_BETA\tGROUP1_SE\tGROUP2_SE\tDIFF_SE\tGROUP1_P\tGROUP2_P\tDIFF_P\tSTD_DIFF_Z\tSTD_DIFF_P\n";
close $fh;
my $config_dir="$dir/new/nested configs";
my $spec={source_mode=>'precomputed_diff_stdized',input_stdized=>$source,
 reference_build=>'hg19',workdir=>"$Bin/..",output_dir=>"$dir/results",
 configs_dir=>$config_dir,artifact_stem=>'path_test',project_tag=>'PATH_TEST',
 gtf_cache_dir=>"$dir/custom gtf cache",top_hit_ld_source=>'HAPLOREG4',
 groups=>[{tag=>'F',files=>[]},{tag=>'M',files=>[]}],
 pairs=>[{pair_tag=>'F_vs_M',group1=>'F',group2=>'M',prefix=>'SEX',label=>'Sex'}]};
open my $sf,'>',"$dir/spec.json" or die $!;
print {$sf} encode_json($spec);close $sf;
my $rc=system($^X,"$Bin/../auto_prepare_and_run_diff_gwas.pl",'--spec',"$dir/spec.json",'--skip-plots','--list-steps');
is($rc,0,'precomputed config generation succeeds in a new nested directory with spaces');
if ($rc==0) {
 my $preset=read_json("$config_dir/auto_path_test_preset.json");
 my $runner=read_json("$config_dir/auto_path_test_runner.json");
 is($preset->{input},$source,'extractor uses existing precomputed input rather than a nonexistent output');
 is($runner->{SOURCE_LONG_GZ},$source,'runner uses existing precomputed input');
 is($runner->{GTF_CACHE_DIR},$spec->{gtf_cache_dir},'SAS runner receives explicit GTF cache');
}
done_testing();
sub read_json {
 open my $fh,'<',$_[0] or die $!;
 my $value=decode_json(do {local $/;<$fh>});close $fh;return $value;
}
