/* Public GTEx v8 tissue medians, not donor-level expression measurements.
   Required original macro definitions are loaded by the reader catalog. */
data work.cg_tpm;
  length gene $16;
  infile "expression40.tsv" dlm='09'x dsd firstobs=2 truncover;
  input gene :$16. Lung Blood Cortex Liver Muscle Heart Pancreas Skin Kidney Spleen;
run;

data work.cg_log2 work.cg_rowz;
  set work.cg_tpm;
  array tissue[10] Lung Blood Cortex Liver Muscle Heart Pancreas Skin Kidney Spleen;
  do i=1 to dim(tissue);
    tissue[i]=log2(tissue[i]+1);
  end;
  output work.cg_log2;
  row_mean=mean(of tissue[*]); row_sd=std(of tissue[*]);
  if n(of tissue[*]) ne 10 or row_sd<=0 then do;
    put 'ERROR: Incomplete or constant gene profile: ' gene=;
    stop;
  end;
  do i=1 to dim(tissue);
    tissue[i]=(tissue[i]-row_mean)/row_sd;
  end;
  output work.cg_rowz;
  keep gene Lung Blood Cortex Liver Muscle Heart Pancreas Skin Kidney Spleen;
run;

ods html5(id=saspy_internal) options(bitmap_mode="inline") device=png;
ods listing gpath="&BOOK_OUT" image_dpi=400 style=journal;
title;footnote;
%global rnd;
%clustergram4sas(dsdin=work.cg_rowz,rowname_var=gene,
  numeric_vars=Lung Blood Cortex Liver Muscle Heart Pancreas Skin Kidney Spleen,
  switch_rownames_as_colnames=1,cluster_type=3,stdize=0,
  perform_log2_for_numeric_vars=0,
  height=1050,width=950,columnweights=0.17 0.83,rowweights=0.14 0.86,
  rowaxis_rotation=vertical,
  colaxis_font_setting=%str(style=normal size=12 family="Albany AMT"),
  rowaxis_font_setting=%str(style=normal size=12 family="Albany AMT"),
  show_heatmap_grid=0,outputfmt=png,heatmap_dsd=work.cg_plot,
  colormodel=,
  rangemap_setting=%str(rangeattrmap name='cgscale';
    range -3 - 3 / rangecolormodel=(CX2166AC CXF7F7F7 CXB2182B);
    endrangeattrmap;
    rangeattrvar var=dist attrmap='cgscale' attrvar=RangeVar;),
  continuouslegend_setting=%str(location=outside valign=bottom border=false
    title="Within-gene Z of log2(median TPM + 1)" valuecounthint=7));

data work.cg_tree;set work.xaxis_dendrogram;run;
data work.cg_tissue_tree;set work.yaxis_dendrogram;run;
proc template;source HeatDendrogram / file="&BOOK_OUT/cluster40_template_original.sas";run;
/* Keep all 40 gene names. The macro's default axis policy thins these ticks.
   Change the generated display template only, not the original macro or data. */
data _null_;
  infile "&BOOK_OUT/cluster40_template_original.sas" lrecl=32767 truncover;
  file "&BOOK_OUT/cluster40_template.sas" lrecl=32767;
  input;
  _infile_=tranwrd(_infile_,'displaysecondary=(tickvalues)',
    'displaysecondary=(tickvalues) discreteopts=(tickvaluefitpolicy=none)');
  put _infile_;
run;
proc template;
  %include "&BOOK_OUT/cluster40_template.sas";
run;
ods listing gpath="&BOOK_OUT" image_dpi=400 style=journal;
ods graphics / reset=index imagename='cluster40' imagefmt=png width=6.4in height=7.1in;
proc sgrender data=work.cg_plot template=HeatDendrogram;run;
ods graphics / reset=index imagename='cluster40' imagefmt=svg width=6.4in height=7.1in;
proc sgrender data=work.cg_plot template=HeatDendrogram;run;

/* Label the existing gene tree; no new clustering is performed. */
ods graphics / reset=index imagename='gene_branch_review' imagefmt=png width=10in height=4.8in;
/* This helper writes no value when it finds no parents to remove. Initialize
   its shared scratch name explicitly in this new caller. */
%global removed_parents;
%let removed_parents=;
%review_tree_branches(inputdsd=cg_tree,
  y_name_var=x_name_,y_parent_var=x_parent_,y_height_var=x_height_,
  outdsd=cg_review,branch_name_dsd=cg_order);
proc template;source HeatDendrogram1 / file="&BOOK_OUT/branch_review_template.sas";run;
/* Suppress duplicate parent text at the zero-height leaves. Higher branch
   labels still use the original macro's generated coordinates and heights. */
data work.cg_review_display;
  set work.cg_review;
  if x_height_=0 then all_parents='';
run;
ods graphics / reset=index imagename='gene_branch_review_clean' imagefmt=png width=10in height=4.8in;
proc sgrender data=work.cg_review_display template=HeatDendrogram1;run;
ods graphics / reset=index imagename='gene_branch_review' imagefmt=svg width=10in height=4.8in;
proc sgrender data=work.cg_review_display template=HeatDendrogram1;run;

/* The original review macro materializes leaf membership for each branch. */
proc sql;
  create table work.cg_branches as
  select all_parents as branch length=16,
    count(distinct x_name_) as genes
  from work.ab_res
  where not prxmatch('/^CL[0-9]+$/',strip(x_name_))
  group by all_parents
  order by genes descending,branch;
quit;
