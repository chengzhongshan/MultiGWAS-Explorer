/*
Run this in SAS ODA after uploading:
  a wide differential GWAS subset .tsv with the expected beta / SE / P columns

This version avoids FILENAME GZIP because the current SAS ODA session rejected
the GZIP filename option.
*/

options mprint mlogic symbolgen;

%include "~/Manhattan4DiffGWASs_png.sas";

__WIDE_IMPORT_BLOCK__

proc sort data=scz_mh;
  by CHR BP;
run;

/* Bind GPLOT symbol colors to chromosome values even if some chromosomes
   are absent from this genome-wide subset. */
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
  fig_width=1200,
  fig_height=__MANHATTAN_FIG_HEIGHT__,
  fontsize=2.6,
  y_axis_label_size=1.8,
  y_axis_value_size=1.8,
  gwas_label_names=%str(__MANHATTAN_GWAS_LABEL_NAMES__),
  gwas_label_x_pct=50,
  gwas_label_y_frac=0.90,
  gwas_label_size=1.8,
  gwas_label_halo_size=1.8,
  gwas_label_angle=0,
  flip1stGWAS_signal=0,
  rm_signals_with_logP_lt=0.5,
  outputfigname=__OUTPUT_PREFIX__,
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
