/*
Run this in SAS ODA after uploading:
  1) a wide differential GWAS subset .tsv.gz with the expected beta / SE / P columns
  2) Manhattan4DiffGWASs_png.sas
 
The input table was prepared by extract_significant_diff_gwas.pl from:
  PGC_SCZ_female_vs_male_diff_effects.stdized.tsv.gz
*/

*options mprint mlogic symbolgen;

ods _all_ close;
ods listing;

%include "~/Manhattan4DiffGWASs_png.sas";

__WIDE_IMPORT_BLOCK__

__MANHATTAN_SORT_BLOCK__

/* GPLOT numbers symbols by the chromosomes actually present. Build that
   symbol list from chromosome values so a missing chromosome does not shift
   the colors of all subsequent chromosomes. */
proc sort data=scz_mh(keep=CHR where=(not missing(CHR)))
          out=_genome_chr_palette nodupkey;
  by CHR;
run;
data _null_;
  set _genome_chr_palette end=last;
  length colors $4096 color $8;
  retain colors '';
  color=scan('cx0072bd cxd95319 cxedb120 cx7e2f8e cx77ac30 cx4dbeee cxa2142f',
             mod(CHR-1,7)+1, ' ');
  colors=catx(' ',colors,color);
  if last then call symputx('genome_chr_colors',colors);
run;

/* Font controls for the Manhattan figure.
   Keep comments outside the macro call so SAS does not misparse the
   keyword-argument list. */
%Manhattan4DiffGWASs(
  dsdin=scz_mh,
  pos_var=BP,
  chr_var=CHR,
  P_var=__MANHATTAN_P_VAR__,
  Other_P_vars=__MANHATTAN_OTHER_P_VARS__,
  logP=1,
  gwas_thrsd=7.30103,
  dotsize=1,
  uniq_colors=&genome_chr_colors,
  _logP_topval=10,
  y_axix_step=2,
  fig_width=__MANHATTAN_FIG_WIDTH__,
  fig_height=__MANHATTAN_FIG_HEIGHT__,
  fontsize=__MANHATTAN_FONTSIZE__,
  y_axis_label_size=__MANHATTAN_Y_AXIS_LABEL_SIZE__,
  y_axis_value_size=__MANHATTAN_Y_AXIS_VALUE_SIZE__,
  gwas_label_names=%str(__MANHATTAN_GWAS_LABEL_NAMES__),
  gwas_label_x_pct=50,
  gwas_label_y_frac=0.90,
  gwas_label_size=__MANHATTAN_GWAS_LABEL_SIZE__,
  gwas_label_halo_size=__MANHATTAN_GWAS_LABEL_HALO_SIZE__,
  gwas_label_angle=0,
  flip1stGWAS_signal=0,
  rm_signals_with_logP_lt=0.5,
  outputfigname=__OUTPUT_PREFIX__,
  plot_title=Genome-wide Manhattan: -log10(P) by chromosome color,
  Use_scaled_pos=1,
  sep_chr_grp=0,
  gwas_sortedby_numchrpos=1
);

data _null_;
  file "~/__OUTPUT_PREFIX___png.html" lrecl=32767;
  put '<!doctype html>';
  put '<html><head><meta charset="utf-8">';
  put '<title>__HTML_TITLE__</title>';
  put '<style>body{margin:0;padding:16px;font-family:Arial,sans-serif;background:#fff;} img{max-width:100%;height:auto;display:block;}</style>';
  put '</head><body>';
  put '<img src="__OUTPUT_PREFIX__.png" alt="__HTML_TITLE__">';
  put '</body></html>';
run;

ods listing close;
