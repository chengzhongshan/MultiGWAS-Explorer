/*
Run this in SAS ODA after uploading:
  1) a requested top-hit CSV with wide beta / SE / P columns
  2) RandBetween.sas
  3) mkfmt4grps_by_var.sas
  4) beta2OR_forest_plot.sas
*/

ods _all_ close;
ods listing;

%include "__RAND_INCLUDE_PATH__";
%include "__MKFMT_INCLUDE_PATH__";
%include "__FOREST_MACRO_INCLUDE_PATH__";

%let forest_track_ids=__FOREST_TRACK_IDS__;
%let forest_track_labels=__FOREST_TRACK_LABELS__;
%let forest_track_beta_vars=__FOREST_TRACK_BETA_VARS__;
%let forest_track_se_vars=__FOREST_TRACK_SE_VARS__;
%let forest_track_p_vars=__FOREST_TRACK_P_VARS__;
%let forest_track_count=__FOREST_TRACK_COUNT__;
%let forest_fig_width=__FOREST_FIG_WIDTH__;
%let forest_fig_height=__FOREST_FIG_HEIGHT__;
%let forest_dotsize=__FOREST_DOTSIZE__;
%let forest_y_font_size=__FOREST_Y_FONT_SIZE__;
%let forest_min_axis=__FOREST_MIN_AXIS__;
%let forest_max_axis=__FOREST_MAX_AXIS__;
%let forest_xaxis_value_range=__FOREST_XAXIS_VALUE_RANGE__;
%let forest_default_hit_class=__FOREST_DEFAULT_HIT_CLASS__;

proc import datafile="__INPUT_CSV_PATH__"
  dbms=csv
  out=forest_hits_raw
  replace;
  getnames=yes;
  guessingrows=max;
run;

data forest_hits;
  set forest_hits_raw;
  length forest_snp $128 forest_gene $256 forest_hit_class $32;
  forest_snp=strip(coalescec(SNP, ''));
  forest_gene=strip(coalescec(gene, 'NA'));
  if missing(forest_gene) then forest_gene='NA';
  forest_hit_class=upcase(strip(coalescec(hit_class, "&forest_default_hit_class")));
  if missing(hit_order) then hit_order=_n_;
  if missing(forest_snp) then delete;
run;

proc datasets lib=work nolist;
  delete rendered_panels;
quit;

data rendered_panels;
  length track_order 8 track_id $64 track_label $128 png_file $256;
  stop;
run;

proc sql noprint;
  select count(*) into: forest_total_hits trimmed
  from forest_hits;
quit;

data _null_;
  length forest_snp_local $128 forest_gene_local $256 forest_title_local $384;
  if symget('forest_total_hits') = '' then call symputx('forest_total_hits', 0, 'L');
  call symputx('forest_single_snp', '', 'L');
  call symputx('forest_single_gene', '', 'L');
  call symputx('forest_single_title', '', 'L');
  set forest_hits(obs=1);
  forest_snp_local = strip(forest_snp);
  forest_gene_local = strip(forest_gene);
  forest_title_local = forest_snp_local;
  if not missing(forest_gene_local) and upcase(forest_gene_local) ne 'NA' then
    forest_title_local = cats(forest_snp_local, ' (', forest_gene_local, ')');
  call symputx('forest_single_snp', forest_snp_local, 'L');
  call symputx('forest_single_gene', forest_gene_local, 'L');
  call symputx('forest_single_title', forest_title_local, 'L');
run;

%macro render_single_snp_forest;
  data forest_single_panel;
    set forest_hits(obs=1);
    length cohort_id $64 cohort_label $128 forest_gene_label $256;
    forest_gene_label=forest_gene;
    %do __single_i=1 %to &forest_track_count;
      cohort_id="%scan(%superq(forest_track_ids), &__single_i, |)";
      cohort_label="%scan(%superq(forest_track_labels), &__single_i, |)";
      forest_order=&__single_i;
      beta_value=input(strip(vvalue(%scan(%superq(forest_track_beta_vars), &__single_i, |))), best32.);
      se_value=input(strip(vvalue(%scan(%superq(forest_track_se_vars), &__single_i, |))), best32.);
      p_value=input(strip(vvalue(%scan(%superq(forest_track_p_vars), &__single_i, |))), best32.);
      if not missing(beta_value) and not missing(se_value) and not missing(p_value) then output;
    %end;
    keep cohort_id cohort_label forest_gene_label forest_order beta_value se_value p_value;
  run;

  proc sql noprint;
    select count(*) into: forest_single_n trimmed
    from forest_single_panel;
  quit;

  %if %sysevalf(%superq(forest_single_n)=,boolean) %then %let forest_single_n=0;
  %if &forest_single_n<=0 %then %return;

  %Beta2OR_forest_plot(
    dsdin=forest_single_panel,
    beta_var=beta_value,
    se_var=se_value,
    sig_p_var=p_value,
    marker_var=cohort_label,
    marker_label=Cohort,
    svgoutname=__OUTPUT_IMAGE_PREFIX_PATH___single_snp,
    plot_title=%superq(forest_single_title),
    figfmt=png,
    figwidth=&forest_fig_width,
    figheight=&forest_fig_height,
    dotsize=&forest_dotsize,
    autolegend=0,
    sort_var4y=forest_order,
    both_y_font_size=&forest_y_font_size,
    sig_datalabel_pos=center,
    sig_datalabel_size=%sysevalf(&forest_y_font_size+1),
    min_axis=&forest_min_axis,
    max_axis=&forest_max_axis,
    xaxis_value_range=%str(&forest_xaxis_value_range),
    randomize_output_suffix=0,
    extra_condition4updatedsd=%nrstr(
      length sigtag $10.;
      grp=1;
      if &sig_p_var<5e-8 and &sig_p_var>0 then sigtag='*';
      else sigtag='';
    ),
    outdsd=forest_single_panel_out
  );

  data _panel_meta;
    length track_order 8 track_id $64 track_label $128 png_file $256;
    track_order=1;
    track_id="single_snp";
    track_label="&forest_single_title";
    png_file="__OUTPUT_IMAGE_BASENAME_PREFIX___single_snp.png";
    output;
  run;

  proc append base=rendered_panels data=_panel_meta force;
  run;
%mend;

%macro render_forest_panel(track_order=, track_id=, track_label=, beta_var=, se_var=, p_var=);
  data forest_panel_pre;
    set forest_hits;
    length forest_marker $128 forest_gene_label $256 forest_hit_class2 $32;
    forest_marker=forest_snp;
    forest_gene_label=forest_gene;
    forest_hit_class2=forest_hit_class;
    forest_hit_order=input(strip(vvalue(hit_order)), best32.);
    if missing(forest_hit_order) then forest_hit_order=_n_;
    beta_value=input(strip(vvalue(&beta_var)), best32.);
    se_value=input(strip(vvalue(&se_var)), best32.);
    p_value=input(strip(vvalue(&p_var)), best32.);
    if missing(beta_value) or missing(se_value) or missing(p_value) then delete;
    if forest_hit_class2='DIFFERENTIAL' then class_order=1;
    else if forest_hit_class2='COMMON' then class_order=2;
    else class_order=3;
    keep forest_marker forest_gene_label forest_hit_class2 forest_hit_order
         beta_value se_value p_value class_order;
  run;

  proc sort data=forest_panel_pre;
    by class_order forest_hit_order forest_marker;
  run;

  data forest_panel;
    set forest_panel_pre;
    by class_order;
    retain forest_order 0;
    forest_order+1;
    if first.class_order and forest_order>1 then separator_value=forest_order-0.5;
  run;

  proc sql noprint;
    select count(*) into: forest_panel_n trimmed
    from forest_panel;
    select put(separator_value, best12.-L) into: forest_sep_values separated by ' '
    from forest_panel
    where not missing(separator_value);
  quit;

  %if %sysevalf(%superq(forest_panel_n)=,boolean) %then %let forest_panel_n=0;
  %if &forest_panel_n<=0 %then %return;

  %Beta2OR_forest_plot(
    dsdin=forest_panel,
    beta_var=beta_value,
    se_var=se_value,
    sig_p_var=p_value,
    marker_var=forest_marker,
    marker_label=SNP,
    svgoutname=__OUTPUT_IMAGE_PREFIX_PATH___&track_id,
    plot_title=&track_label,
    figfmt=png,
    figwidth=&forest_fig_width,
    figheight=&forest_fig_height,
    dotsize=&forest_dotsize,
    autolegend=0,
    sort_var4y=forest_order,
    y2axis_ticket_var=forest_gene_label,
    both_y_font_size=&forest_y_font_size,
    sig_datalabel_pos=center,
    sig_datalabel_size=%sysevalf(&forest_y_font_size+1),
    min_axis=&forest_min_axis,
    max_axis=&forest_max_axis,
    xaxis_value_range=%str(&forest_xaxis_value_range),
    y_refline_values=%superq(forest_sep_values),
    randomize_output_suffix=0,
    extra_condition4updatedsd=%nrstr(
      length sigtag $10.;
      if forest_hit_class2='COMMON' then grp=1;
      else if forest_hit_class2='DIFFERENTIAL' then grp=2;
      else grp=3;
      if &sig_p_var<5e-8 and &sig_p_var>0 then do;
        sigtag='*';
        grp=grp+10;
      end;
      else sigtag='';
    ),
    outdsd=forest_panel_out
  );

  data _panel_meta;
    length track_order 8 track_id $64 track_label $128 png_file $256;
    track_order=&track_order;
    track_id="&track_id";
    track_label="&track_label";
    png_file=cats("__OUTPUT_IMAGE_BASENAME_PREFIX___", strip("&track_id"), ".png");
    output;
  run;

  proc append base=rendered_panels data=_panel_meta force;
  run;
%mend;

%macro render_all_forest_panels;
  %local i track_id track_label beta_var se_var p_var;
  %do i=1 %to &forest_track_count;
    %let track_id=%qscan(%superq(forest_track_ids), &i, |);
    %let track_label=%qscan(%superq(forest_track_labels), &i, |);
    %let beta_var=%scan(%superq(forest_track_beta_vars), &i, |);
    %let se_var=%scan(%superq(forest_track_se_vars), &i, |);
    %let p_var=%scan(%superq(forest_track_p_vars), &i, |);
    %render_forest_panel(
      track_order=&i,
      track_id=&track_id,
      track_label=&track_label,
      beta_var=&beta_var,
      se_var=&se_var,
      p_var=&p_var
    );
  %end;
%mend;

data _null_;
  length _render_cmd $64;
  _forest_total_hits = input(symget('forest_total_hits'), best32.);
  if missing(_forest_total_hits) then _forest_total_hits = 0;
  if _forest_total_hits = 1 then _render_cmd = '%nrstr(%render_single_snp_forest;)';
  else _render_cmd = '%nrstr(%render_all_forest_panels;)';
  call execute(_render_cmd);
run;

proc sort data=rendered_panels;
  by track_order track_id;
run;

data _null_;
  file "__OUTPUT_HTML_PATH__" lrecl=32767;
  put '<!doctype html>';
  put '<html><head><meta charset="utf-8">';
  put '<title>__FOREST_HTML_TITLE__</title>';
  put '<style>';
  put 'body{margin:0;padding:20px;font-family:Arial,sans-serif;background:#fff;color:#1f2937;}';
  put '.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:24px;align-items:start;}';
  put '.panel{margin:0;padding:14px;border:1px solid #d9dee7;border-radius:10px;background:#fff;box-shadow:0 1px 4px rgba(0,0,0,0.06);}';
  put '.panel figcaption{font-size:18px;font-weight:700;margin-bottom:10px;text-align:center;}';
  put '.panel img{width:100%;height:auto;display:block;}';
  put '</style></head><body>';
  put '<div class="grid">';
  do until(eof);
    set rendered_panels end=eof;
    put '<figure class="panel">';
    put '<figcaption>' track_label +(-1) '</figcaption>';
    put '<img src="' png_file +(-1) '" alt="' track_label +(-1) ' forest plot">';
    put '</figure>';
  end;
  put '</div></body></html>';
run;

data _null_;
  file "__OUTPUT_MANIFEST_PATH__" lrecl=32767;
  put 'track_order' '09'x 'track_id' '09'x 'track_label' '09'x 'png_file';
  do until(eof);
    set rendered_panels end=eof;
    put track_order best12. '09'x track_id '09'x track_label '09'x png_file;
  end;
run;

ods listing close;
