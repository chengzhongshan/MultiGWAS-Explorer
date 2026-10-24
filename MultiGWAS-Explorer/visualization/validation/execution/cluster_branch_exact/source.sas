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

/* Branch selection follows the same full-tree execution. */
/* Continue after cluster40.sas. CL10 is an observed branch of this exact
   40-gene input and source version, not a stable biological identifier. */
%let chosen_branch=CL10;
proc sql noprint;
  select count(*) into :branch_exists trimmed
  from work.cg_tree where x_name_="&chosen_branch";
quit;
%macro require_branch;
  %if &branch_exists ne 1 %then %do;
    %put ERROR: Requested branch does not exist uniquely.;
    %abort cancel;
  %end;
%mend;
%require_branch;

ods listing gpath="&BOOK_OUT" image_dpi=300 style=journal;
ods graphics / reset=index imagename='selected_tree' imagefmt=png width=6.4in height=3.5in;
%subbranches(data=cg_tree,
  y_name_var=x_name_,y_parent_var=x_parent_,y_height_var=x_height_,
  out=work.cg_subtree,parent=&chosen_branch,height=,render=YES);

/* Preserve the native result. Its reconstructed root passes through a macro
   variable and can lose numeric precision. Reattach original numeric heights
   by OLD_X_NAME_, without changing node membership or topology. */
data work.cg_subtree_native;set work.cg_subtree;run;
proc sql;
  create table work.cg_subtree_exact as
  select a.*,b.x_height_ as source_height
  from work.cg_subtree as a left join work.cg_tree as b
    on a.old_x_name_=b.x_name_;
quit;
data work.cg_subtree;
  set work.cg_subtree_exact;
  x_height_=source_height;
  drop source_height;
run;

proc sql;
  create table work.cg_members as
  select distinct x_name_ as gene length=16
  from work.cg_subtree
  where not prxmatch('/^CL[0-9]+$/',strip(x_name_));
  create table work.cg_selected_values as
  select a.* from work.cg_rowz as a
  inner join work.cg_members as b on a.gene=b.gene
  order by a.gene;
quit;

/* Verify the selected tree keeps its original node heights. SUBBRANCHES
   renumbers CL nodes and retains their old names in OLD_X_NAME_. */
proc sql;
  create table work.cg_selected_heights as
  select old_x_name_ as original_node length=32,x_height_ as height
  from work.cg_subtree order by original_node;
  create table work.cg_expected_heights as
  select a.x_name_ as original_node length=32,a.x_height_ as height
  from work.cg_tree as a
  inner join work.cg_selected_heights as b
    on a.x_name_=b.original_node
  order by original_node;
quit;
/* A diagnostic comparison records the native root discrepancy separately. */
proc sql;
  create table work.cg_native_heights as
  select old_x_name_ as original_node length=32,x_height_ as height
  from work.cg_subtree_native order by original_node;
quit;
proc compare base=work.cg_expected_heights compare=work.cg_native_heights
  criterion=1e-12 method=absolute noprint;
  id original_node;var height;
run;
%let native_compare_sysinfo=&sysinfo;
proc compare base=work.cg_expected_heights compare=work.cg_selected_heights
  criterion=1e-12 method=absolute noprint;
  id original_node;var height;
run;
%let branch_compare_sysinfo=&sysinfo;
data work.cg_branch_check;
  length comparison $40;
  comparison='Selected node heights unchanged';
  sysinfo=&branch_compare_sysinfo;
  native_sysinfo=&native_compare_sysinfo;
run;

/* Reuse the generated clustergram template with the selected ORIGINAL tree.
   No PROC CLUSTER call is made on this eight-gene subset. */
proc sql;
  create table work.cg_branch_cells as
  select a.row,a.col,a.dist
  from work.cg_plot as a inner join work.cg_members as b
    on a.col=b.gene;
quit;
data work.cg_x_selected;
  set work.cg_subtree;keep x_name_ x_parent_ x_height_;
run;
data work.cg_branch_plot;
  merge work.cg_branch_cells work.cg_x_selected work.cg_tissue_tree;
run;
ods listing gpath="&BOOK_OUT" image_dpi=400 style=journal;
ods graphics / reset=index imagename='cluster_selected' imagefmt=png width=6.4in height=4.2in;
proc sgrender data=work.cg_branch_plot template=HeatDendrogram;run;
ods graphics / reset=index imagename='cluster_selected' imagefmt=svg width=6.4in height=4.2in;
proc sgrender data=work.cg_branch_plot template=HeatDendrogram;run;
