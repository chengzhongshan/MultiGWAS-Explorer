#!/usr/bin/env bash
# Cygwin usage:
# bash examples/run_real_gtex_boxplots.sh [MATRIX_OR_TWO_GENE_TSV] [ODA_METADATA_DIR]
# The input may be the full local matrix or an already extracted two-gene TSV.
# Only the small selected expression table is uploaded. Real GTEx UMAP metadata
# must already exist online; no full expression matrix is downloaded by SAS.
set -euo pipefail
export PATH=/usr/bin:$PATH
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
default_matrix='E:/LongCOVID_HGI_GWAS/GTEx_Exp_SAS_Data/exp.matrix/exprMatrix_ensgSymbolCombined.tsv'
if command -v cygpath >/dev/null 2>&1; then
  default_matrix="$(cygpath -u "$default_matrix")"
else
  default_matrix='/mnt/e/LongCOVID_HGI_GWAS/GTEx_Exp_SAS_Data/exp.matrix/exprMatrix_ensgSymbolCombined.tsv'
fi
matrix="${1:-$default_matrix}"
metadata_dir="${2:-~/data/GTEx_V8}"
[[ -f "$matrix" ]] || { printf 'Missing local expression matrix: %s\n' "$matrix" >&2; exit 1; }
matrix="$(cd "$(dirname "$matrix")" && pwd -P)/$(basename "$matrix")"
cd "$project_dir"
output_dir="$(mktemp -d "$project_dir/SAS_Output_gtex_cache_real_XXXXXX")"
printf 'Output directory: %s\n' "$output_dir"
perl tools/subset_gtex_expression.pl --input "$matrix" \
  --genes 'MAP3K19 R3HDM1' --output "$output_dir/gtex_two_genes.tsv"
upload="$output_dir/gtex_real_${output_dir##*_}.tsv"

# Convert just two wide gene rows into typed sample/value records. Preserve the
# raw numeric values exactly; SAS applies log2(value+1) once, after the join.
row_count="$(perl - "$output_dir/gtex_two_genes.tsv" "$upload" "$output_dir/input_manifest.json" "$matrix" <<'PERL'
use strict; use warnings; use JSON::PP; use Digest::SHA qw(sha256_hex);
my ($input,$output,$manifest,$source)=@ARGV;
open my $in,'<:raw',$input or die "$input: $!\n";
my $header=<$in>; die "Empty input\n" unless defined $header;
$header=~s/\r?\n\z//;
my @samples=split /\t/,$header,-1; shift @samples;
my %seen;
for (@samples) {die "Invalid or repeated sample ID\n" unless /^GTEX-[A-Za-z0-9-]+$/ && !$seen{$_}++;}
open my $out,'>:raw',$output or die "$output: $!\n";
print {$out} "rownames\tsample_id\traw_exp\n";
my (%genes,$rows); $rows=0;
while (my $line=<$in>) {
  $line=~s/\r?\n\z//; my ($gene,@values)=split /\t/,$line,-1;
  my $symbol=(split /\|/,$gene)[-1];
  die "Unexpected or repeated gene\n" unless $symbol=~/^(?:MAP3K19|R3HDM1)$/ && !$genes{$symbol}++;
  die "Incorrect sample count\n" unless @values==@samples;
  for my $i (0..$#samples) {
    my $v=$values[$i];
    die "Invalid expression for $gene $samples[$i]\n" unless $v=~/^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
    print {$out} "$gene\t$samples[$i]\t$v\n"; ++$rows;
  }
}
close $in or die $!; close $out or die $!;
die "Both genes are required\n" unless keys(%genes)==2;
open my $hash,'<:raw',$input or die $!;
my $sha=Digest::SHA->new(256)->addfile($hash)->hexdigest; close $hash;
open my $audit,'>:raw',$manifest or die $!;
print {$audit} JSON::PP->new->canonical->pretty->encode({source_matrix=>$source,
  subset_sha256=>$sha,sample_columns=>scalar(@samples),expression_rows=>$rows,
  genes=>[sort keys %genes],transform=>'SAS: log2(raw_exp+1)',
  phenotype_join=>'sample_id = UMAP.seq_id; unmatched/duplicate IDs fail',
  outlier_rule=>'exclude clusters containing OUTLIER, matching the original macro'});
close $audit or die $!; print $rows;
PERL
)"
remote_input="~/${upload##*/}"
metadata_sas="${metadata_dir//\'/\'\'}"
{
  cat examples/sas_macros/Boxplots4GenesInGTExV8_Cached.sas
  printf '\n%%let real_input=%s;\n%%let real_expected_n=%s;\n%%let real_metadata_dir=%s;\n' \
    "$remote_input" "$row_count" "$metadata_sas"
  cat examples/sas/gtex_real_prepare_and_plot.sas
} > "$output_dir/real_gtex_boxplot.sas"
printf '\nSubmitting %s expression rows. Exact command:\n' "$row_count"
printf 'perl run_sas_codes_or_file_in_ODA.pl --upload-file "%s" --file "%s" --output-prefix "%s" --run-timeout-seconds 300\n' \
  "$upload" "$output_dir/real_gtex_boxplot.sas" "$output_dir/plot"
OPEN_RESULT=0 SAS_ODA_AUTOLOAD_MACROS=0 \
  perl run_sas_codes_or_file_in_ODA.pl --upload-file "$upload" \
    --file "$output_dir/real_gtex_boxplot.sas" --output-prefix "$output_dir/plot" \
    --run-timeout-seconds 300 2>&1 | tee "$output_dir/command.log"
printf '\nReal-data validation:\n'
grep -E '^GTEX_REAL_PASS:' "$output_dir/plot/output.html.info.txt"
printf '\nHTML results:\n'
find "$output_dir/plot" -type f -name '*.html' -print
printf '\nFor a faster rerun, use the saved two-gene subset:\n'
printf 'bash examples/run_real_gtex_boxplots.sh "%s" "%s"\n' "$output_dir/gtex_two_genes.tsv" "$metadata_dir"
