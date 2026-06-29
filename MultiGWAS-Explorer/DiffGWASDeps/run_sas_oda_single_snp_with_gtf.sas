/*
Run this in SAS ODA after uploading:
  1) a wide differential GWAS subset .tsv.gz with the expected beta / SE / P columns

This script draws a local Manhattan plot with gene tracks for a single target SNP by:
  - reading the uploaded wide GWAS subset
  - computing Z-scores from beta / SE columns
  - locating the requested target SNP
  - subsetting a configurable +/- local window around that SNP
  - ensuring the requested GTF dataset exists in SAS ODA
  - calling SNP_Local_Manhattan_With_GTF with configurable track settings
*/

*options mprint mlogic symbolgen;

%let target_snp=__TARGET_SNP__;
%let local_window_bp=__LOCAL_WINDOW_BP__;
%let gtf_label_snps=__GTF_LABEL_SNPS__;
%let html_outfile=__OUTPUT_HTML__;
%let gwas_dsd=__GWAS_DATASET__;
%let target_hit_dsd=__TARGET_HIT_DATASET__;
%let target_local_dsd=__TARGET_LOCAL_DATASET__;

%let gtf_dsd=__GTF_DSD__;
%let fm_libpath=__FM_LIBPATH__;
%let gtf_local_dsd=__GTF_LOCAL_DSD__;
%let gtf_gz_url=__GTF_GZ_URL__;
%let gtf_assoc_pvars=__GTF_ASSOC_PVARS__;
%let gtf_zscore_vars=__GTF_ZSCORE_VARS__;
%let gtf_labels=__GTF_LABELS__;
%let gtf_dist2snp=__GTF_DIST2SNP__;
%let gtf_design_width=__GTF_DESIGN_WIDTH__;
%let gtf_design_height=__GTF_DESIGN_HEIGHT__;
%let gtf_dist2sep_genes=__GTF_DIST2SEP_GENES__;
%let gtf_shift_text_yval=__GTF_SHIFT_TEXT_YVAL__;
%let gtf_pct4neg_y=__GTF_PCT4NEG_Y__;
%let gtf_adjval4header=__GTF_ADJVAL4HEADER__;
%let gtf_yaxis_label=__GTF_YAXIS_LABEL__;
%let gtf_colorbar_label=__GTF_COLORBAR_LABEL__;
%let gtf_yoffset4textlabels=__GTF_YOFFSET4TEXTLABELS__;
%let gtf_include_non_protein_coding=__GTF_INCLUDE_NON_PROTEIN_CODING__;

ods _all_ close;
ods html5 file="~/&html_outfile"
  options(bitmap_mode='inline')
  style=HTMLBlue;
ods graphics on / outputfmt=png;

%if %length(&fm_libpath) > 0 %then %do;
  libname FM "&fm_libpath";
%end;

%include "~/adj_grpnum4close_gene_bed_regs.sas";
%include "~/Multgscatter_with_gene_exons.sas";
%include "~/map_grp_assoc2gene4covidsexgwas.sas";
%include "~/SNP_Local_Manhattan_With_GTF.sas";
%include "~/Lattice_gscatter_over_bed_track.sas";

__GTF_IMPORT_BLOCK__

%macro _ensure_effective_gtf_dsd(target_chr=,target_bp=,flank_bp=);
  %global effective_gtf_dsd effective_gtf_dist2snp effective_gtf_pct4neg_y effective_gtf_dist2sep_genes effective_gtf_grp_font_size;
  %local _n_genes_for_plot _n_exons_for_plot _n_non_protein_for_plot _n_tx_for_plot _compact_gtf_dsd;
  %let effective_gtf_dsd=&gtf_local_dsd;
  %let effective_gtf_dist2snp=&gtf_dist2snp;
  %let effective_gtf_pct4neg_y=&gtf_pct4neg_y;
  %let effective_gtf_dist2sep_genes=&gtf_dist2sep_genes;
  %let effective_gtf_grp_font_size=8;
  %if %sysevalf(&local_window_bp > &effective_gtf_dist2snp) %then %do;
    %put NOTE: Expanding local GTF half-window from &effective_gtf_dist2snp bp to the requested SNP-centered local window (&local_window_bp bp).;
    %let effective_gtf_dist2snp=&local_window_bp;
  %end;

  %let _compact_gtf_dsd=&effective_gtf_dsd._compact;

  proc sql;
    create table _gtf_transcript_rank as
    select coalescec(genesymbol,gene) as gene_key length=256,
           coalescec(transcript_id,transcript_name) as tx_key length=256,
           sum(type='exon') as exon_count,
           (max(en)-min(st)+1) as tx_span
    from &effective_gtf_dsd
    where type in ('transcript','exon')
      and not missing(coalescec(genesymbol,gene))
      and not missing(coalescec(transcript_id,transcript_name))
    group by calculated gene_key, calculated tx_key
    ;
  quit;

  proc sort data=_gtf_transcript_rank;
    by gene_key descending exon_count descending tx_span tx_key;
  run;

  data _gtf_transcript_best;
    set _gtf_transcript_rank;
    by gene_key;
    if first.gene_key;
  run;

  proc sql;
    create table &_compact_gtf_dsd as
    select a.*
    from &effective_gtf_dsd as a
    left join _gtf_transcript_best as b
      on coalescec(a.genesymbol,a.gene)=b.gene_key
    where a.type='gene'
       or (
            a.type in ('transcript','exon')
        and not missing(coalescec(a.transcript_id,a.transcript_name))
        and coalescec(a.transcript_id,a.transcript_name)=b.tx_key
       )
    ;
  quit;

  %if %sysfunc(exist(&_compact_gtf_dsd)) %then %let effective_gtf_dsd=&_compact_gtf_dsd;

  proc sql noprint;
    select count(distinct genesymbol),
           sum(type='exon'),
           sum(original_protein_coding=0),
           count(distinct coalescec(transcript_id,transcript_name))
      into :_n_genes_for_plot trimmed,
           :_n_exons_for_plot trimmed,
           :_n_non_protein_for_plot trimmed,
           :_n_tx_for_plot trimmed
    from &effective_gtf_dsd
    ;
  quit;

  %put NOTE: Using local GTF dataset &effective_gtf_dsd for single-SNP plotting.;
  %put NOTE: Local GTF stats: genes=&_n_genes_for_plot representative_transcripts=&_n_tx_for_plot exons=&_n_exons_for_plot non_protein_coding_features=&_n_non_protein_for_plot include_non_protein_coding=&gtf_include_non_protein_coding.;
%mend;

%macro _auto_tune_gene_track_ratio(target_chr=,target_bp=,gtf_dsd=);
  %global effective_gtf_pct4neg_y effective_gtf_dist2sep_genes effective_gtf_grp_font_size;
  %local _nearby_gene_count _recommended_pct4neg_y _recommended_dist2sep_genes _recommended_grp_font_size _display_window_bp _normalized_base_dist2sep_genes;
  %let _nearby_gene_count=0;
  %let _recommended_pct4neg_y=&gtf_pct4neg_y;
  %let _normalized_base_dist2sep_genes=&gtf_dist2sep_genes;
  %let _recommended_dist2sep_genes=&gtf_dist2sep_genes;
  %let _recommended_grp_font_size=&effective_gtf_grp_font_size;
  %let _display_window_bp=%sysevalf(2*&effective_gtf_dist2snp);
  %if %sysevalf(&_normalized_base_dist2sep_genes >= %sysevalf(&_display_window_bp*0.25)) %then %do;
    %let _normalized_base_dist2sep_genes=%sysfunc(max(50000,%sysfunc(int(%sysevalf(&_display_window_bp*0.005)))));
  %end;
  %let _recommended_dist2sep_genes=&_normalized_base_dist2sep_genes;

  proc sql noprint;
    select count(distinct coalescec(genesymbol,gene))
      into :_nearby_gene_count trimmed
    from &gtf_dsd
    where not missing(strip(cats(gene)))
      and type='gene'
      and chr=&target_chr
      and st <= (&target_bp + &effective_gtf_dist2snp)
      and en >= (&target_bp - &effective_gtf_dist2snp)
    ;
  quit;

  %if %sysevalf(%superq(_nearby_gene_count)=,boolean) %then %let _nearby_gene_count=0;

  %if %sysevalf(&_nearby_gene_count > 120) %then %do;
    %let _recommended_pct4neg_y=%sysfunc(max(&gtf_pct4neg_y,3.00));
    %let _recommended_dist2sep_genes=%sysfunc(min(&_normalized_base_dist2sep_genes,%sysfunc(max(50000,%sysfunc(int(%sysevalf(&_display_window_bp*0.002)))))));
    %let _recommended_grp_font_size=5;
  %end;
  %else %if %sysevalf(&_nearby_gene_count > 80) %then %do;
    %let _recommended_pct4neg_y=%sysfunc(max(&gtf_pct4neg_y,2.40));
    %let _recommended_dist2sep_genes=%sysfunc(min(&_normalized_base_dist2sep_genes,%sysfunc(max(75000,%sysfunc(int(%sysevalf(&_display_window_bp*0.003)))))));
    %let _recommended_grp_font_size=5;
  %end;
  %else %if %sysevalf(&_nearby_gene_count > 50) %then %do;
    %let _recommended_pct4neg_y=%sysfunc(max(&gtf_pct4neg_y,1.90));
    %let _recommended_dist2sep_genes=%sysfunc(min(&_normalized_base_dist2sep_genes,%sysfunc(max(100000,%sysfunc(int(%sysevalf(&_display_window_bp*0.004)))))));
    %let _recommended_grp_font_size=6;
  %end;
  %else %if %sysevalf(&_nearby_gene_count > 30) %then %do;
    %let _recommended_pct4neg_y=%sysfunc(max(&gtf_pct4neg_y,1.50));
    %let _recommended_dist2sep_genes=%sysfunc(min(&_normalized_base_dist2sep_genes,%sysfunc(max(125000,%sysfunc(int(%sysevalf(&_display_window_bp*0.006)))))));
    %let _recommended_grp_font_size=6;
  %end;
  %else %if %sysevalf(&_nearby_gene_count > 15) %then %do;
    %let _recommended_pct4neg_y=%sysfunc(max(&gtf_pct4neg_y,1.20));
    %let _recommended_dist2sep_genes=%sysfunc(min(&_normalized_base_dist2sep_genes,%sysfunc(max(150000,%sysfunc(int(%sysevalf(&_display_window_bp*0.008)))))));
    %let _recommended_grp_font_size=7;
  %end;

  %let effective_gtf_pct4neg_y=&_recommended_pct4neg_y;
  %let effective_gtf_dist2sep_genes=&_recommended_dist2sep_genes;
  %let effective_gtf_grp_font_size=&_recommended_grp_font_size;
  %put NOTE: Auto-tuned single-SNP gene-track ratio using nearby_gene_count=&_nearby_gene_count pct4neg_y=&effective_gtf_pct4neg_y (base=&gtf_pct4neg_y) dist2sep_genes=&effective_gtf_dist2sep_genes (base=&gtf_dist2sep_genes) grp_font_size=&effective_gtf_grp_font_size.;
%mend;

%macro _maybe_expand_gtf_dist2snp(target_chr=,target_bp=,gtf_dsd=);
  %local _nearby_gene_count;
  %let _nearby_gene_count=0;
  proc sql noprint;
    select count(*) into: _nearby_gene_count trimmed
    from &gtf_dsd
    where not missing(strip(cats(gene)))
      and type='gene'
      and chr=&target_chr
      and (&target_bp between (st-&gtf_dist2snp) and (en+&gtf_dist2snp))
    ;
  quit;

  %if %sysevalf(&_nearby_gene_count=0) %then %do;
    %let effective_gtf_dist2snp=&local_window_bp;
    %put NOTE: No genes were found within the default GTF distance of &gtf_dist2snp bp. Expanding gene-track search to the local window size (&effective_gtf_dist2snp bp).;
  %end;
%mend;

__WIDE_IMPORT_BLOCK__

proc sort data=&gwas_dsd;
  by CHR BP;
run;

proc sql;
  create table &target_hit_dsd as
  select *
  from &gwas_dsd
  where strip(SNP)=strip("&target_snp");
quit;

proc sql outobs=1 noprint;
  select CHR, BP into :target_chr trimmed, :target_bp trimmed
  from &target_hit_dsd;
quit;

%if %superq(target_chr)= or %superq(target_bp)= %then %do;
  %put ERROR: Target SNP &target_snp was not found in the uploaded GWAS subset.;
  %abort 255;
%end;

%_ensure_effective_gtf_dsd(target_chr=&target_chr,target_bp=&target_bp,flank_bp=&local_window_bp);
%_maybe_expand_gtf_dist2snp(target_chr=&target_chr,target_bp=&target_bp,gtf_dsd=&effective_gtf_dsd);
%_auto_tune_gene_track_ratio(target_chr=&target_chr,target_bp=&target_bp,gtf_dsd=&effective_gtf_dsd);

proc sql;
  create table &target_local_dsd as
  select *
  from &gwas_dsd
  where CHR=&target_chr
    and BP between (&target_bp-&local_window_bp) and (&target_bp+&local_window_bp)
  order by CHR, BP;
quit;

%SNP_Local_Manhattan_With_GTF(
  gwas_dsd=&target_local_dsd,
  chr_var=CHR,
  AssocPVars=&gtf_assoc_pvars,
  SNP_IDs=&target_snp,
  dist2snp=&effective_gtf_dist2snp,
  SNP_Var=SNP,
  Pos_Var=BP,
  gtf_dsd=&effective_gtf_dsd,
  ZscoreVars=&gtf_zscore_vars,
  gwas_labels_in_order=&gtf_labels,
  design_width=&gtf_design_width,
  design_height=&gtf_design_height,
  barthickness=10,
  dotsize=5,
  grp_font_size=&effective_gtf_grp_font_size,
  dist2sep_genes=&effective_gtf_dist2sep_genes,
  where_cndtn_for_gwasdsd=%str(),
  shift_text_yval=&gtf_shift_text_yval,
  fig_fmt=png,
  pct4neg_y=&effective_gtf_pct4neg_y,
  adjval4header=&gtf_adjval4header,
  scatter_yaxis_label=&gtf_yaxis_label,
  heatmap_legend_title=&gtf_colorbar_label,
  makedotheatmap=1,
  makeheatmapdotintooneline=0,
  SNPs2label_scatterplot_dots=&gtf_label_snps,
  Yoffset4textlabels=&gtf_yoffset4textlabels,
  verbose=0
);

ods html5 close;
