/* Local desktop-SAS debug paths emitted automatically by the pipeline. */
%let local_debug_root=\cygdrive\c\Users\zcheng\Downloads\AI4Coding\perlMCP4AOA_GWAS\MultiGWAS-Explorer\MultiGWAS-Explorer;
%let local_debug_deps=&local_debug_root\DiffGWASDeps;
%let local_debug_output_dir=&local_debug_root;
%let local_debug_output_html=\cygdrive\c\Users\zcheng\Downloads\AI4Coding\perlMCP4AOA_GWAS\MultiGWAS-Explorer\MultiGWAS-Explorer\AOA_GWAS_Data_diff_SAS_manhattan_png.html;
%let local_debug_manhattan_macro=\cygdrive\c\Users\zcheng\Downloads\AI4Coding\perlMCP4AOA_GWAS\MultiGWAS-Explorer\MultiGWAS-Explorer\Manhattan4DiffGWASs_png.local_debug.sas;
%let local_debug_wide_data_gz=\cygdrive\c\Users\zcheng\Downloads\AI4Coding\perlMCP4AOA_GWAS\MultiGWAS-Explorer\MultiGWAS-Explorer\cache\sas_manhattan\AOA_GWAS_Data_diff_SAS_manhattan.p_lt_0_05.tsv.gz;


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

%include "&local_debug_manhattan_macro";

filename mhdata zip "&local_debug_wide_data_gz" gzip;

data scz_mh;
  infile mhdata dlm='09'x dsd firstobs=2 truncover lrecl=32767;
  input
    CHR
    BP
    MP2PRT_DS_ALL_DIFF_P
    DS_ALL_P
    MP2PRT_P
    META_P
  ;
run;


/* Compact Manhattan input was sorted by numeric CHR and BP locally. */


/* Font controls for the Manhattan figure.
   Keep comments outside the macro call so SAS does not misparse the
   keyword-argument list. */
%Manhattan4DiffGWASs(
  dsdin=scz_mh,
  pos_var=BP,
  chr_var=CHR,
  P_var=MP2PRT_DS_ALL_DIFF_P,
  Other_P_vars=DS_ALL_P MP2PRT_P META_P,
  logP=1,
  gwas_thrsd=7.30103,
  dotsize=1,
  _logP_topval=10,
  y_axix_step=2,
  fig_width=1800,
  fig_height=700,
  fontsize=2.4,
  y_axis_label_size=2.4,
  y_axis_value_size=2.2,
  gwas_label_names=%str(MP2PRT vs DS_ALL raw differential P|DS_ALL association P|MP2PRT association P|Meta association P),
  gwas_label_x_pct=50,
  gwas_label_y_frac=0.90,
  gwas_label_size=2.4,
  gwas_label_halo_size=2.4,
  gwas_label_angle=0,
  flip1stGWAS_signal=0,
  rm_signals_with_logP_lt=0.5,
  outputfigname=AOA_GWAS_Data_diff_SAS_manhattan,
  Use_scaled_pos=1,
  sep_chr_grp=0,
  gwas_sortedby_numchrpos=1
);

data _null_;
  file "&local_debug_output_html" lrecl=32767;
  put '<!doctype html>';
  put '<html><head><meta charset="utf-8">';
  put '<title>AOA_GWAS_Data_diff differential GWAS Manhattan Plot</title>';
  put '<style>body{margin:0;padding:16px;font-family:Arial,sans-serif;background:#fff;} img{max-width:100%;height:auto;display:block;}</style>';
  put '</head><body>';
  put '<img src="AOA_GWAS_Data_diff_SAS_manhattan.png" alt="AOA_GWAS_Data_diff differential GWAS Manhattan Plot">';
  put '</body></html>';
run;

ods listing close;
