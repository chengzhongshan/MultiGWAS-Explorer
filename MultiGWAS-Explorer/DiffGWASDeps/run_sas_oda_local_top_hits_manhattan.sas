/*
Run this in SAS ODA after uploading:
  1) a wide differential GWAS subset .tsv.gz with the expected beta / SE / P columns
  2) Manhattan4DiffGWASs_png.sas

This script extends the genome-wide differential GWAS workflow by selecting
top loci, annotating them with nearby genes from HaploReg, subsetting +/- window
regions around those loci, and drawing the stacked local Manhattan plot only.

For the gene-track local plot that uses SNP_Local_Manhattan_With_GTF, use the
companion runner/script pair:
  - run_sas_oda_local_top_hits_with_gtf.sas
  - run_sas_oda_local_top_hits_with_gtf_download_html.sh
*/

*options mprint mlogic symbolgen;

%let top_hit_focus_pvar=__TOP_HIT_FOCUS_PVAR__;
%let top_hit_mode=__TOP_HIT_MODE__;
%let top_hit_filter_expr=__TOP_HIT_FILTER_EXPR__;
%let top_hit_signal_thrshd=__TOP_HIT_SIGNAL_THRSHD__;
%let top_hit_signal_thrshds=__TOP_HIT_SIGNAL_THRSHDS__;
%let top_hit_dist_bp=__TOP_HIT_DIST_BP__;
%let top_hit_selection_method=__TOP_HIT_SELECTION_METHOD__;
%let top_hit_ld_r2_threshold=__TOP_HIT_LD_R2_THRESHOLD__;
%let top_hit_ld_populations=__TOP_HIT_LD_POPULATIONS__;
%let top_hit_ld_population_rule=__TOP_HIT_LD_POPULATION_RULE__;
%let top_hit_ld_query_failure_action=__TOP_HIT_LD_QUERY_FAILURE_ACTION__;
%let top_hit_ld_cache_basename=__TOP_HIT_LD_CACHE_BASENAME__;
%let top_hit_ld_cache_min_r2=__TOP_HIT_LD_CACHE_MIN_R2__;
%let top_hit_ld_cache_miss_action=__TOP_HIT_LD_CACHE_MISS_ACTION__;
%let top_hit_max_loci=__TOP_HIT_MAX_LOCI__;
%let top_hit_ld_audit_basename=__TOP_HIT_LD_AUDIT_BASENAME__;
%let local_window_bp=__LOCAL_WINDOW_BP__;
%let local_max_hits_per_fig=__LOCAL_MAX_HITS_PER_FIG__;
%let local_n_gwas_tracks=%eval(1 + %sysfunc(countw(%str(__MANHATTAN_OTHER_P_VARS__),%str( ))));
%let local_top_hits_csv_basename=__LOCAL_TOP_HITS_CSV_BASENAME__;
%let requested_top_hits_csv_basename=__REQUESTED_TOP_HITS_CSV_BASENAME__;
%let target_snp_list=__TARGET_SNP_LIST__;
%let target_snp_gene_map=__TARGET_SNP_GENES__;
%let local_ld_snps=__GTF_LD_SNPS__;
%let local_ld_marker_symbol=__GTF_LD_MARKER_SYMBOL__;
%let local_ld_marker_color=__GTF_LD_MARKER_COLOR__;
%let local_ld_display_mode=__GTF_LD_DISPLAY_MODE__;
%let local_ld_r2_values=__GTF_LD_R2_VALUES__;
%let local_ld_heatmap_colors=__GTF_LD_HEATMAP_COLORS__;
%let local_ld_heatmap_legend_title=__GTF_LD_HEATMAP_LEGEND_TITLE__;
%let common_assoc_pvars=__COMMON_ASSOC_P_VARS__;
%let lmh_angle4xaxis_label=__LOCAL_MANHATTAN_ANGLE4XAXIS_LABEL__;
%let lmh_xgrp_y_pos=__LOCAL_MANHATTAN_XGRP_Y_POS__;
%let lmh_yoffset_top=__LOCAL_MANHATTAN_YOFFSET_TOP__;
%let lmh_yoffset_bottom=__LOCAL_MANHATTAN_YOFFSET_BOTTOM__;
%let lmh_fontsize=__LOCAL_MANHATTAN_FONTSIZE__;
%let lmh_y_axis_label_size=__LOCAL_MANHATTAN_Y_AXIS_LABEL_SIZE__;
%let lmh_y_axis_value_size=__LOCAL_MANHATTAN_Y_AXIS_VALUE_SIZE__;
%let gtf_dsd=__GTF_DSD__;
%let fm_libpath=__FM_LIBPATH__;
%let gtf_local_dsd=__GTF_LOCAL_DSD__;
%let gtf_gz_url=__GTF_GZ_URL__;
%let gtf_include_non_protein_coding=__GTF_INCLUDE_NON_PROTEIN_CODING__;

ods _all_ close;
ods listing;

%macro _initialize_optional_library;
%if %length(&fm_libpath) > 0 %then %do;
  libname FM "&fm_libpath";
%end;
%mend;
%_initialize_optional_library;

%include "~/get_top_signal_within_dist.sas";
%include "~/QueryLD_SNPs_at_Haploreg4.sas";
%include "~/get_top_signal_with_ld.sas";
%include "~/__GET_GTF_MACRO_BASENAME__";
%include "~/Manhattan4DiffGWASs_png.sas";

%macro _load_optional_ld_cache;
%if not %sysevalf(%superq(top_hit_ld_cache_basename)=,boolean) %then %do;
proc import datafile="~/&top_hit_ld_cache_basename"
  out=top_hit_ld_cache dbms=tab replace;
  guessingrows=max;
run;
proc datasets library=work nolist;
  modify top_hit_ld_cache;
  index create query_population=(query_snp ld_population);
quit;
%let HAPLOREG_LD_CACHE_DSD=work.top_hit_ld_cache;
%let HAPLOREG_LD_CACHE_MIN_R2=&top_hit_ld_cache_min_r2;
%let HAPLOREG_LD_CACHE_MISS_ACTION=&top_hit_ld_cache_miss_action;
%put NOTE: Loaded candidate-specific HaploReg LD cache &top_hit_ld_cache_basename (minimum r2=&top_hit_ld_cache_min_r2).;
%end;
%mend;
%_load_optional_ld_cache;

%macro _find_first_column(lib=,mem=,outvar=,candidates=);
  %global &outvar;
  %local _i _cand _matched;
  %let &outvar=;
  %do _i=1 %to %sysfunc(countw(%superq(candidates),%str( )));
    %let _cand=%scan(%superq(candidates),&_i,%str( ));
    %let _matched=;
    proc sql noprint;
      select name into: _matched trimmed
      from dictionary.columns
      where libname="%upcase(&lib)"
        and memname="%upcase(&mem)"
        and upcase(name)="%upcase(&_cand)"
      ;
    quit;
    %if %superq(_matched) ne %then %do;
      %let &outvar=&_matched;
      %goto _done_find_first_column;
    %end;
  %end;
  %_done_find_first_column:
%mend;

%macro _load_requested_target_snps(outdsd=);
  %global requested_target_snps_loaded;
  %let requested_target_snps_loaded=0;

  %if %sysevalf(%superq(target_snp_list)=,boolean) %then %return;

  data &outdsd;
    length SNP $128 _raw $32767;
    _raw=prxchange('s/[[:space:]]+/,/o', -1, strip(symget('target_snp_list')));
    _raw=prxchange('s/,+/,/o', -1, _raw);
    _raw=prxchange('s/^,+|,+$//o', -1, _raw);
    do hit_order=1 to countw(_raw, ',', 'm');
      SNP=strip(scan(_raw, hit_order, ',', 'm'));
      if not missing(SNP) then output;
    end;
    keep hit_order SNP;
  run;

  proc sort data=&outdsd nodupkey;
    by hit_order SNP;
  run;

  proc sql noprint;
    select count(*) into: requested_target_snps_loaded trimmed
    from &outdsd
    ;
  quit;

  %if %sysevalf(%superq(requested_target_snps_loaded)=,boolean) %then %let requested_target_snps_loaded=0;
  %put NOTE: Loaded &requested_target_snps_loaded requested target SNP(s) from TARGET_SNP_LIST.;
%mend;

%macro _load_requested_target_snp_genes(outdsd=);
  %global req_target_snp_genes_loaded;
  %let req_target_snp_genes_loaded=0;

  %if %sysevalf(%superq(target_snp_gene_map)=,boolean) %then %return;

  data &outdsd;
    length SNP $128 gene $256 _raw _entry $32767;
    _raw=prxchange('s/[[:space:]]+/,/o', -1, strip(symget('target_snp_gene_map')));
    _raw=prxchange('s/,+/,/o', -1, _raw);
    _raw=prxchange('s/^,+|,+$//o', -1, _raw);
    do _i=1 to countw(_raw, ',', 'm');
      _entry=strip(scan(_raw, _i, ',', 'm'));
      if not missing(_entry) then do;
        SNP=strip(scan(_entry, 1, ':'));
        gene=strip(substr(_entry, lengthn(SNP) + 2));
        if not missing(SNP) and not missing(gene) then output;
      end;
    end;
    keep SNP gene;
  run;

  proc sort data=&outdsd nodupkey;
    by SNP;
  run;

  proc sql noprint;
    select count(*) into: req_target_snp_genes_loaded trimmed
    from &outdsd
    ;
  quit;

  %if %sysevalf(%superq(req_target_snp_genes_loaded)=,boolean) %then %let req_target_snp_genes_loaded=0;
  %put NOTE: Loaded &req_target_snp_genes_loaded requested target SNP gene override(s) from TARGET_SNP_GENE_MAP.;
%mend;

%macro _load_req_top_hits_csv(outdsd=);
  %global requested_top_hits_loaded;
  %let requested_top_hits_loaded=0;

  %if %superq(requested_top_hits_csv_basename)= %then %return;

  proc import datafile="~/&requested_top_hits_csv_basename"
    out=&outdsd
    dbms=csv
    replace;
    guessingrows=max;
  run;

  %if not %sysfunc(exist(work.&outdsd)) %then %return;

  %_find_first_column(lib=WORK,mem=%upcase(&outdsd),outvar=req_chr_var,candidates=CHR chr);
  %_find_first_column(lib=WORK,mem=%upcase(&outdsd),outvar=req_bp_var,candidates=BP bp POS pos position);
  %_find_first_column(lib=WORK,mem=%upcase(&outdsd),outvar=req_snp_var,candidates=SNP rsid marker id);
  %_find_first_column(lib=WORK,mem=%upcase(&outdsd),outvar=req_hit_order_var,candidates=hit_order HIT_ORDER order panel_order rank);
  %_find_first_column(lib=WORK,mem=%upcase(&outdsd),outvar=req_gene_var,candidates=gene GENE nearest_gene gene_symbol symbol);
  %_find_first_column(lib=WORK,mem=%upcase(&outdsd),outvar=req_snp_gene_var,candidates=snp_gene SNP_GENE snp2gene);

  %if %superq(req_chr_var)= or %superq(req_bp_var)= or %superq(req_snp_var)= %then %do;
    %put WARNING: Requested top-hit CSV lacks resolvable CHR/BP/SNP columns, so it will be ignored.;
    %put WARNING: Resolved columns: chr=&req_chr_var bp=&req_bp_var snp=&req_snp_var hit_order=&req_hit_order_var;
    proc datasets library=work nolist;
      delete &outdsd;
    quit;
    %return;
  %end;

  data &outdsd;
    set &outdsd;
    length CHR 8 BP 8 SNP $128 hit_order 8 gene $256 snp_gene $128;
    CHR=input(strip(vvaluex("&req_chr_var")),best32.);
    BP=input(strip(vvaluex("&req_bp_var")),best32.);
    SNP=strip(vvaluex("&req_snp_var"));
    %if %superq(req_hit_order_var) ne %then %do;
    hit_order=input(strip(vvaluex("&req_hit_order_var")),best32.);
    %end;
    %else %do;
    hit_order=_n_;
    %end;
    %if %superq(req_gene_var) ne %then %do;
    gene=strip(vvaluex("&req_gene_var"));
    %end;
    %else %do;
    gene='';
    %end;
    %if %superq(req_snp_gene_var) ne %then %do;
    snp_gene=strip(vvaluex("&req_snp_gene_var"));
    %end;
    %else %do;
    snp_gene='';
    %end;
    if upcase(gene)='NA' then gene='';
    if upcase(snp_gene)='NA' then snp_gene='';
    if missing(CHR) or missing(BP) or missing(SNP) then delete;
  run;

  proc sort data=&outdsd nodupkey;
    by hit_order CHR BP SNP;
  run;

  proc sql noprint;
    select count(*) into: requested_top_hits_loaded trimmed
    from &outdsd
    ;
  quit;

  %if %sysevalf(%superq(requested_top_hits_loaded)=,boolean) %then %let requested_top_hits_loaded=0;
%mend;

%macro _prepare_gtf_gene_fallback(top_hits_dsd=,gtf_dsd=,outdsd=);
  %local _gtf_lib _gtf_mem;
  %let _gtf_lib=%upcase(%scan(%superq(gtf_dsd),1,.));
  %let _gtf_mem=%upcase(%scan(%superq(gtf_dsd),2,.));
  %if %superq(_gtf_mem)= %then %do;
    %let _gtf_mem=&_gtf_lib;
    %let _gtf_lib=WORK;
  %end;

  %_find_first_column(lib=&_gtf_lib,mem=&_gtf_mem,outvar=gtf_chr_var,candidates=chr seqname chromosome chrom chr_raw);
  %_find_first_column(lib=&_gtf_lib,mem=&_gtf_mem,outvar=gtf_start_var,candidates=start st bp1 txstart);
  %_find_first_column(lib=&_gtf_lib,mem=&_gtf_mem,outvar=gtf_end_var,candidates=end en bp2 txend);
  %_find_first_column(lib=&_gtf_lib,mem=&_gtf_mem,outvar=gtf_gene_var,candidates=gene gene_name gene_symbol symbol name2 name gene_id transcript_name transcript_id);
  %_find_first_column(lib=&_gtf_lib,mem=&_gtf_mem,outvar=gtf_feature_var,candidates=feature type);

  %if %superq(gtf_chr_var)= or %superq(gtf_start_var)= or %superq(gtf_end_var)= or %superq(gtf_gene_var)= %then %do;
    %put WARNING: Could not resolve enough GTF columns in &gtf_dsd for fallback gene mapping.;
    %put WARNING: Resolved columns: chr=&gtf_chr_var start=&gtf_start_var end=&gtf_end_var gene=&gtf_gene_var feature=&gtf_feature_var;
    data &outdsd;
      length rsid $40 gtf_gene $256;
      stop;
    run;
    %return;
  %end;

  proc sql;
    create table _gtf_gene_fallback_candidates as
    select a.SNP as rsid length=40,
           b.&gtf_gene_var as gtf_gene length=256,
           case
             when a.BP between b.&gtf_start_var and b.&gtf_end_var then 0
             when a.BP < b.&gtf_start_var then (b.&gtf_start_var - a.BP)
             else (a.BP - b.&gtf_end_var)
           end as gtf_gene_distance
    from &top_hits_dsd as a, &gtf_dsd as b
    where prxmatch('/^rs[0-9]+$/i', strip(a.SNP)) > 0
      and not missing(strip(cats(b.&gtf_gene_var)))
      and b.&gtf_start_var <= b.&gtf_end_var
      and upcase(compress(strip(cats(a.CHR)),'CHR')) = upcase(compress(strip(cats(b.&gtf_chr_var)),'CHR'))
      %if %superq(gtf_feature_var) ne %then %do;
      and (missing(strip(cats(b.&gtf_feature_var))) or lowcase(strip(cats(b.&gtf_feature_var)))='gene')
      %end;
    ;
  quit;

  proc sort data=_gtf_gene_fallback_candidates;
    by rsid gtf_gene_distance gtf_gene;
  run;

  data &outdsd;
    set _gtf_gene_fallback_candidates;
    by rsid;
    if first.rsid;
    keep rsid gtf_gene;
  run;
%mend;

%macro _set_gtf_region_lists(top_hits_dsd=,flank_bp=);
  %global gtf_region_chrs gtf_region_starts gtf_region_ends;
  %let gtf_region_chrs=;
  %let gtf_region_starts=;
  %let gtf_region_ends=;
  proc sql noprint;
    select strip(put(CHR,best32.)),
           strip(put(case when BP>&flank_bp then BP-&flank_bp else 1 end,best32.)),
           strip(put(BP+&flank_bp,best32.))
      into :gtf_region_chrs separated by '|',
           :gtf_region_starts separated by '|',
           :gtf_region_ends separated by '|'
    from &top_hits_dsd
    where prxmatch('/^rs[0-9]+$/i', strip(SNP)) > 0
    order by CHR, BP
    ;
  quit;
%mend;

%macro _ensure_effective_gtf_dsd(top_hits_dsd=);
  %global effective_gtf_dsd;
  %local _n_gtf_genes _n_gtf_non_protein;
  %let effective_gtf_dsd=gtf_local_mh_plot;

  %_set_gtf_region_lists(top_hits_dsd=&top_hits_dsd,flank_bp=&local_window_bp);
  %if %superq(gtf_region_chrs)= %then %do;
    %put WARNING: No rsID-style top-hit regions were available for region-limited GTF extraction.;
    %let effective_gtf_dsd=&gtf_dsd;
    %return;
  %end;

  %__GET_GTF_MACRO_NAME__(
    gtf_gz_url=&gtf_gz_url,
    outdsd=&gtf_local_dsd,
    region_chrs=&gtf_region_chrs,
    region_starts=&gtf_region_starts,
    region_ends=&gtf_region_ends
  );

  data &effective_gtf_dsd;
    length chr 8 chr_text $64 ensembl $64 type $32 genesymbol gene $256 protein_coding 8 original_protein_coding 8 _bio_text $32767;
    set &gtf_local_dsd(rename=(chr=chr_text));
    chr_raw=coalescec(chr_raw,seqname,chr_text);
    chr_text=upcase(prxchange('s/^CHR//i',1,strip(coalescec(chr_text,chr_raw,seqname))));
    if chr_text='X' then chr=23;
    else if chr_text='Y' then chr=24;
    else if chr_text in ('M','MT') then chr=25;
    else chr=input(chr_text,best32.);
    if missing(chr) then delete;
    if missing(ensembl) then ensembl=coalescec(source,'gencode');
    type=lowcase(coalescec(type,feature));
    if missing(st) then st=start;
    if missing(en) then en=end;
    gene=coalescec(gene,gene_name,gene_id,transcript_name,transcript_id,feature);
    genesymbol=coalescec(gene_name,gene,transcript_name,gene_id,transcript_id,feature);
    _bio_text=lowcase(coalescec(gene_type,transcript_type,attribute,''));
    original_protein_coding=(index(_bio_text,'protein_coding')>0);
    if missing(protein_coding) then protein_coding=original_protein_coding;
    %if %sysevalf(&gtf_include_non_protein_coding,boolean) %then %do;
    protein_coding=1;
    %end;
    %else %do;
    protein_coding=original_protein_coding;
    %end;
    if type not in ('gene','transcript','exon') then delete;
    drop _bio_text;
  run;

  proc sql noprint;
    select count(distinct genesymbol),
           sum(original_protein_coding=0)
      into :_n_gtf_genes trimmed,
           :_n_gtf_non_protein trimmed
    from &effective_gtf_dsd
    ;
  quit;

  %put NOTE: Using locally generated GTF dataset &effective_gtf_dsd for fallback top-hit gene labeling.;
  %put NOTE: Local GTF stats: genes=&_n_gtf_genes non_protein_coding_features=&_n_gtf_non_protein include_non_protein_coding=&gtf_include_non_protein_coding.;
%mend;

__WIDE_IMPORT_BLOCK__

proc format library=work;
  value $__bootstrap_fmt__
    0='0'
  ;
run;

proc sort data=scz_mh;
  by CHR BP;
run;

%macro _prepare_req_top_hit_inputs;
%global requested_top_hits_loaded;
%let requested_top_hits_loaded=0;
%if %upcase(&top_hit_mode)=COMMON_ASSOCIATION %then %do;
data scz_mh;
  set scz_mh;
  COMMON_ASSOC_P=min(of &common_assoc_pvars);
run;
%end;

%_load_requested_target_snps(outdsd=requested_target_snps);
%if %sysevalf(&requested_target_snps_loaded,boolean) %then %do;
%let requested_top_hits_loaded=0;
data requested_top_hits_csv;
  length CHR 8 BP 8 SNP $128 hit_order 8 gene $256 snp_gene $128;
  stop;
run;
%end;
%else %do;
%_load_req_top_hits_csv(outdsd=requested_top_hits_csv);
%if not %sysfunc(exist(work.requested_top_hits_csv)) %then %do;
data requested_top_hits_csv;
  length CHR 8 BP 8 SNP $128 hit_order 8 gene $256 snp_gene $128;
  stop;
run;
%end;

%if %sysevalf(&requested_top_hits_loaded,boolean)
    and %upcase(&top_hit_mode)=COMMON_ASSOCIATION %then %do;
data requested_top_hits_csv;
  set requested_top_hits_csv;
  COMMON_ASSOC_P=min(of &common_assoc_pvars);
run;
%end;
%end;
%mend;
%_prepare_req_top_hit_inputs;

%macro _pick_top_hits_by_thr(candidate_dsd=);
  %local _i _n _thr _n_hits_this;
  %let _n=%sysfunc(countw(%superq(top_hit_signal_thrshds),%str( )));
  %if &_n=0 %then %let _n=1;
  %do _i=1 %to &_n;
    %let _thr=%scan(%superq(top_hit_signal_thrshds),&_i,%str( ));
    %if %superq(_thr)= %then %let _thr=&top_hit_signal_thrshd;
    %put NOTE: Trying top-hit threshold &_thr for mode=&top_hit_mode focus=&top_hit_focus_pvar;
    %let _n_hits_this=0;
    proc sql noprint;
      select count(*) into: _n_hits_this trimmed
      from &candidate_dsd
      where (&top_hit_focus_pvar>0) and (&top_hit_focus_pvar<&_thr);
    quit;
    %if %sysevalf(&_n_hits_this>0) %then %do;
      %let top_hit_signal_thrshd=&_thr;
      %goto _picked;
    %end;
  %end;
  %_picked:
%mend;

%macro _select_computed_top_hits(candidate_dsd=,outdsd=);
  %if %upcase(&top_hit_selection_method)=LD %then %do;
    %put NOTE: Selecting LD-clumped top hits (r2 >= &top_hit_ld_r2_threshold, populations=&top_hit_ld_populations, rule=&top_hit_ld_population_rule).;
    %get_top_signal_with_ld(
      dsdin=&candidate_dsd,
      snp_var=SNP,
      grp_var=CHR,
      signal_var=&top_hit_focus_pvar,
      select_smallest_signal=1,
      pos_var=BP,
      signal_thrshd=&top_hit_signal_thrshd,
      ld_pops=&top_hit_ld_populations,
      ld_r2_threshold=&top_hit_ld_r2_threshold,
      ld_population_rule=&top_hit_ld_population_rule,
      query_failure_action=&top_hit_ld_query_failure_action,
      fallback_dist_bp=&top_hit_dist_bp,
      max_leads=&top_hit_max_loci,
      dsdout=&outdsd,
      audit_out=top_hit_ld_audit
    );
    data &outdsd;
      set &outdsd;
      requested_hit_order=LD_LEAD_RANK;
    run;
  %end;
  %else %do;
    %put NOTE: Selecting top hits by legacy physical-distance pruning because TOP_HIT_SELECTION_METHOD=&top_hit_selection_method.;
    %get_top_signal_within_dist(
      dsdin=&candidate_dsd,
      grp_var=CHR,
      signal_var=&top_hit_focus_pvar,
      select_smallest_signal=1,
      pos_var=BP,
      pos_dist_thrshd=&top_hit_dist_bp,
      dsdout=&outdsd,
      signal_thrshd=&top_hit_signal_thrshd
    );
    data &outdsd;
      set &outdsd;
      requested_hit_order=.;
      length INDEPENDENCE_METHOD $24;
      INDEPENDENCE_METHOD='DISTANCE';
      LD_R2_THRESHOLD=.;
      LD_FALLBACK_BP=&top_hit_dist_bp;
    run;
  %end;
%mend;

%macro _select_top_hit_candidates;
%if %sysevalf(&requested_target_snps_loaded,boolean) %then %do;
proc sql;
  create table top_hit4diffp_raw as
  select a.*,
         b.hit_order as requested_hit_order
  from scz_mh as a
  inner join requested_target_snps as b
    on upcase(strip(a.SNP))=upcase(strip(b.SNP))
  order by b.hit_order, a.CHR, a.BP, a.SNP
  ;
quit;
data top_hit4diffp_raw;
  set top_hit4diffp_raw;
  length INDEPENDENCE_METHOD $24;
  INDEPENDENCE_METHOD='USER_TARGET';
run;
%end;
%else %do;
data top_hit_candidates;
  set scz_mh;
  if &top_hit_filter_expr;
run;

%if %sysevalf(&requested_top_hits_loaded,boolean) %then %do;
  %_pick_top_hits_by_thr(candidate_dsd=requested_top_hits_csv);
  %_select_computed_top_hits(candidate_dsd=requested_top_hits_csv,outdsd=top_hit4diffp_raw);
  data top_hit4diffp_raw;
    length A1 A2 $8 requested_hit_order 8;
    set top_hit4diffp_raw;
    requested_hit_order=coalesce(LD_LEAD_RANK,hit_order,_n_);
    SNP=strip(SNP);
    if missing(A1) then A1=strip(coalescec(EFFECT_ALLELE, ALTERNATIVE_ALLELE));
    if missing(A2) then A2=strip(coalescec(OTHER_ALLELE, REFERENCE_ALLELE));
  run;
  proc sort data=top_hit4diffp_raw;
    by requested_hit_order CHR BP SNP;
  run;
  %put NOTE: The uploaded MAF-filtered CSV supplied LD-clumping candidates, and only LD-independent leads are retained.;
%end;
%else %do;
  %_pick_top_hits_by_thr(candidate_dsd=top_hit_candidates);
  %_select_computed_top_hits(candidate_dsd=top_hit_candidates,outdsd=top_hit4diffp_raw);
%end;
%end;
%mend;
%_select_top_hit_candidates;

%_load_requested_target_snp_genes(outdsd=requested_target_snp_genes);
%macro _ensure_req_target_snp_genes;
%if not %sysfunc(exist(work.requested_target_snp_genes)) %then %do;
  data requested_target_snp_genes;
    length SNP $128 gene $256;
    stop;
  run;
%end;
%mend;
%_ensure_req_target_snp_genes;

proc sql noprint;
  select strip(SNP) into: top_snps separated by ' '
  from top_hit4diffp_raw
  where not missing(SNP)
  order by coalesce(requested_hit_order, 999999999), CHR, BP;
quit;

%macro _check_top_snps_;
%if %superq(top_snps)= %then %do;
  %put ERROR: No top SNPs passed the local top-hit filter and threshold settings.;
  %abort 255;
%end;
%mend;
%_check_top_snps_;

%macro _build_haploreg_gene_map;
  %QueryHaploreg(
    rsids=&top_snps,
    dsdout=snps2genes,
    print_html=0
  );

  %if %sysfunc(exist(WORK.SNPS2GENES)) %then %do;
    data snps2genes_clean;
      length rsid $40 gene $256;
      set snps2genes;
      rsid=strip(rsid);
      gene=strip(gene);
      if prxmatch('/^rs[0-9]+$/i', rsid)=0 then delete;
      if prxmatch('/<[^>]+>/', gene) > 0 then gene='';
      if prxmatch('/could not connect/i', gene) > 0 then gene='';
      if prxmatch('/^detail view for /i', gene) > 0 then gene='';
      gene=prxchange('s/\s*[,;|].*$//',1,gene);
      if upcase(gene)='NA' then gene='';
      if missing(gene) then delete;
    run;

    proc sort data=snps2genes_clean nodupkey;
      by rsid;
    run;
  %end;
  %else %do;
    data snps2genes_clean;
      length rsid $40 gene $256;
      stop;
    run;
  %end;
%mend;

%macro _prepare_top_hit_gene_map;
%if %sysevalf(&requested_target_snps_loaded,boolean) %then %do;
  /* Explicit target-SNP runs use TARGET_SNP_GENE_MAP below.  Avoid making
     their success depend on optional HaploReg or GTF network lookups. */
  data snps2genes_clean;
    length rsid $40 gene $256;
    stop;
  run;

  data snps2genes_gtf_fallback;
    length rsid $40 gtf_gene $256;
    stop;
  run;
%end;
%else %if %sysevalf(&requested_top_hits_loaded,boolean) %then %do;
  data snps2genes_clean;
    length rsid $40 gene $256;
    set requested_top_hits_csv;
    rsid=strip(SNP);
    gene=strip(gene);
    if missing(gene) then gene=scan(strip(snp_gene),2,':');
    if upcase(gene)='NA' then gene='';
    if missing(rsid) then delete;
    keep rsid gene;
  run;

  proc sort data=snps2genes_clean nodupkey;
    by rsid;
  run;

  data snps2genes_gtf_fallback;
    length rsid $40 gtf_gene $256;
    stop;
  run;
%end;
%else %do;
  %_build_haploreg_gene_map;

  %_ensure_effective_gtf_dsd(top_hits_dsd=top_hit4diffp_raw);

  %_prepare_gtf_gene_fallback(
    top_hits_dsd=top_hit4diffp_raw,
    gtf_dsd=&effective_gtf_dsd,
    outdsd=snps2genes_gtf_fallback
  );
%end;
%mend;
%_prepare_top_hit_gene_map;

proc sql;
  create table top_hit4diffp as
  select a.*,
         coalescec(u.gene, r.gene, b.gene, c.gtf_gene) as gene length=256,
         case
           when not missing(u.gene) then 'USER'
           when not missing(r.gene) then 'REQUESTED'
           when not missing(b.gene) then 'HaploReg'
           when not missing(c.gtf_gene) then 'GTF'
           else 'NA'
         end as gene_source length=16,
         coalescec(r.snp_gene, catx(':', a.SNP, coalescec(u.gene, r.gene, b.gene, c.gtf_gene, 'NA'))) as snp_gene length=128
  from top_hit4diffp_raw as a
  left join requested_target_snp_genes as u
    on upcase(strip(a.SNP))=upcase(strip(u.SNP))
  left join requested_top_hits_csv as r
    on a.CHR=r.CHR
   and a.BP=r.BP
   and upcase(strip(a.SNP))=upcase(strip(r.SNP))
  left join snps2genes_clean as b
    on a.SNP=b.rsid
  left join snps2genes_gtf_fallback as c
    on a.SNP=c.rsid
  ;

  create table top_local_signals as
  select a.*, catx(':', b.SNP, coalescec(b.gene,'NA')) as snp_gene length=128
  from scz_mh as a, top_hit4diffp as b
  where a.CHR=b.CHR
    and a.BP between (b.BP-&local_window_bp) and (b.BP+&local_window_bp)
  ;

  select catx(':', SNP, coalescec(gene,'NA')) into: snp_gene_label separated by ' '
  from top_hit4diffp
  order by coalesce(requested_hit_order, 999999999), CHR, BP;
quit;

proc sql;
  create table top_hit_groups as
  select a.snp_gene length=128,
         b.CHR,
         b.BP,
         b.requested_hit_order
  from (
    select distinct snp_gene
    from top_local_signals
  ) as a
  inner join (
    select catx(':', SNP, coalescec(gene,'NA')) as snp_gene length=128,
           CHR,
           BP,
           requested_hit_order
    from top_hit4diffp
  ) as b
    on a.snp_gene=b.snp_gene
  order by coalesce(b.requested_hit_order, 999999999), b.CHR, b.BP
  ;

  select count(*) into: n_top_hits trimmed
  from top_hit_groups;
quit;

data top_hit_groups;
  set top_hit_groups;
  hit_order=_n_;
run;

proc sql noprint;
  select name into: top_hit_export_extra_vars separated by ', '
  from dictionary.columns
  where libname='WORK'
    and memname='TOP_HIT4DIFFP'
    and upcase(name) not in (
      'CHR','BP','SNP','GENE','SNP_GENE','A1','A2','REQUESTED_HIT_ORDER'
    )
  order by varnum
  ;
quit;

proc sql;
  create table top_hit4diffp_export as
  select g.hit_order,
         ceil(g.hit_order / &local_max_hits_per_fig) as panel_index,
         t.CHR,
         t.BP,
         t.SNP,
         t.A1 as EFFECT_ALLELE length=8,
         t.A2 as OTHER_ALLELE length=8,
         t.A2 as REFERENCE_ALLELE length=8,
         t.A1 as ALTERNATIVE_ALLELE length=8,
         t.gene,
         t.snp_gene,
         t.&top_hit_focus_pvar as focus_signal
         %if %superq(top_hit_export_extra_vars) ne %then %do;
         , &top_hit_export_extra_vars
         %end;
  from top_hit4diffp as t
  inner join top_hit_groups as g
    on t.snp_gene=g.snp_gene
  order by g.hit_order
  ;
quit;

proc export data=top_hit4diffp_export
  outfile="~/&local_top_hits_csv_basename"
  dbms=csv
  replace;
run;

%macro _export_top_hit_ld_audit;
%if %upcase(&top_hit_selection_method)=LD
    and %sysfunc(exist(work.top_hit_ld_audit))
    and not %sysevalf(%superq(top_hit_ld_audit_basename)=,boolean) %then %do;
proc export data=top_hit_ld_audit
  outfile="~/&top_hit_ld_audit_basename"
  dbms=tab
  replace;
run;
%end;
%mend;
%_export_top_hit_ld_audit;

data _null_;
  length _lmhpf 8 _nth 8 _nhb 8 _ngt 8;
  _lmhpf=input(symget('local_max_hits_per_fig'), best.);
  if missing(_lmhpf) then _lmhpf=30;
  if _lmhpf<1 then _lmhpf=999999;
  _ngt=input(symget('local_n_gwas_tracks'), best.);
  if missing(_ngt) then _ngt=1;
  if _lmhpf>30 then _lmhpf=30;
  _nth=input(symget('n_top_hits'), best.);
  if missing(_nth) then _nth=0;
  if _nth>0 then _nhb=ceil(_nth/_lmhpf);
  else _nhb=0;
  call symputx('local_max_hits_per_fig', strip(put(_lmhpf, best.)));
  call symputx('n_hit_batches', strip(put(_nhb, best.)));
run;
%put NOTE: local_n_gwas_tracks=&local_n_gwas_tracks;
%put NOTE: local_max_hits_per_fig=&local_max_hits_per_fig;
%put NOTE: n_top_hits=&n_top_hits;
%put NOTE: n_hit_batches=&n_hit_batches;

%macro render_local_hit_batches;
  %local batch_idx start_idx end_idx batch_snp_gene_label batch_output_prefix batch_count batch_dsdin;
  %local local_xgrp_y_pos local_yoffset_setting local_gwas_label_y_frac;
  %local local_angle4xaxis_label local_fontsize local_y_axis_label_size local_y_axis_value_size;
  %do batch_idx=1 %to &n_hit_batches;
    %let start_idx=%eval((&batch_idx-1)*&local_max_hits_per_fig + 1);
    %let end_idx=%sysfunc(min(&n_top_hits,%eval(&batch_idx*&local_max_hits_per_fig)));
    %let batch_count=%eval(&end_idx-&start_idx+1);
    proc sql noprint;
      select snp_gene into: batch_snp_gene_label separated by ' '
      from top_hit_groups
      where hit_order between &start_idx and &end_idx
      order by hit_order;
    quit;
    %put NOTE: batch_idx=&batch_idx start_idx=&start_idx end_idx=&end_idx;
    %put NOTE: batch_count=&batch_count local_n_gwas_tracks=&local_n_gwas_tracks;
    %put NOTE: batch_snp_gene_label=&batch_snp_gene_label;

    %if &batch_idx=1 %then %let batch_output_prefix=__OUTPUT_PREFIX__;
    %else %let batch_output_prefix=__OUTPUT_PREFIX___part&batch_idx;
    %let batch_dsdin=top_local_signals_batch&batch_idx;

    %if &local_n_gwas_tracks>=6 and &batch_count>=6 %then %do;
      %let local_xgrp_y_pos=-6.3;
      %let local_yoffset_setting=%str(offset=(19,0.5));
      %let local_gwas_label_y_frac=0.54;
    %end;
    %else %if &local_n_gwas_tracks>=6 and &batch_count>=4 %then %do;
      %let local_xgrp_y_pos=-5.9;
      %let local_yoffset_setting=%str(offset=(18,0.5));
      %let local_gwas_label_y_frac=0.56;
    %end;
    %else %if &local_n_gwas_tracks>=4 and &batch_count>=6 %then %do;
      %let local_xgrp_y_pos=-5.5;
      %let local_yoffset_setting=%str(offset=(17,0.5));
      %let local_gwas_label_y_frac=0.58;
    %end;
    %else %if &local_n_gwas_tracks>=4 %then %do;
      %let local_xgrp_y_pos=-5.0;
      %let local_yoffset_setting=%str(offset=(16,0.5));
      %let local_gwas_label_y_frac=0.62;
    %end;
    %else %if &batch_count>=6 %then %do;
      %let local_xgrp_y_pos=-4.6;
      %let local_yoffset_setting=%str(offset=(16,0.5));
      %let local_gwas_label_y_frac=0.66;
    %end;
    %else %do;
      %let local_xgrp_y_pos=-4.0;
      %let local_yoffset_setting=%str(offset=(15,0.5));
      %let local_gwas_label_y_frac=0.70;
    %end;
    %put NOTE: local_xgrp_y_pos=&local_xgrp_y_pos;
    %put NOTE: local_yoffset_setting=&local_yoffset_setting;
    %put NOTE: local_gwas_label_y_frac=&local_gwas_label_y_frac;

    %let local_angle4xaxis_label=90;
    *This macro var control the font size of x-axis labels included in the annotation area, i.e., for the snp-gene labels;
    %let local_fontsize=1.8;
    %let local_y_axis_label_size=2.0;
    %let local_y_axis_value_size=1.8;

        %if %superq(lmh_angle4xaxis_label) ne %then %let local_angle4xaxis_label=&lmh_angle4xaxis_label;
        %if %superq(lmh_xgrp_y_pos) ne %then %let local_xgrp_y_pos=&lmh_xgrp_y_pos;
        %if %superq(lmh_yoffset_top) ne %then %do;
          %if %superq(lmh_yoffset_bottom) ne %then
            %let local_yoffset_setting=%str(offset=(&lmh_yoffset_top,&lmh_yoffset_bottom));
          %else
            %let local_yoffset_setting=%str(offset=(&lmh_yoffset_top,0.5));
        %end;
        %if %superq(lmh_fontsize) ne %then %let local_fontsize=&lmh_fontsize;
        %if %superq(lmh_y_axis_label_size) ne %then %let local_y_axis_label_size=&lmh_y_axis_label_size;
        %if %superq(lmh_y_axis_value_size) ne %then %let local_y_axis_value_size=&lmh_y_axis_value_size;

    %put NOTE: local_angle4xaxis_label=&local_angle4xaxis_label;
    %put NOTE: local_fontsize=&local_fontsize;
    %put NOTE: local_y_axis_label_size=&local_y_axis_label_size;
    %put NOTE: local_y_axis_value_size=&local_y_axis_value_size;

    proc sql;
      create table &batch_dsdin as
      select a.*
      from top_local_signals as a
      inner join top_hit_groups as b
        on a.snp_gene=b.snp_gene
      where b.hit_order between &start_idx and &end_idx
      ;
    quit;
    %put NOTE: batch_output_prefix=&batch_output_prefix;
    %put NOTE: batch_dsdin=&batch_dsdin;

    %Manhattan4DiffGWASs(
      dsdin=&batch_dsdin,
      pos_var=BP,
      chr_var=snp_gene,
      P_var=__MANHATTAN_P_VAR__,
      Other_P_vars=__MANHATTAN_OTHER_P_VARS__,
      logP=1,
      gwas_thrsd=7.30103,
      thrsd_line_color=gray,
      dotsize=1,
      _logP_topval=8,
      y_axix_step=2,
      fig_width=__MANHATTAN_FIG_WIDTH__,
      fig_height=__MANHATTAN_FIG_HEIGHT__,
      fontsize=&local_fontsize,
      y_axis_label_size=&local_y_axis_label_size,
      y_axis_value_size=&local_y_axis_value_size,
      gwas_label_names=%str(__MANHATTAN_GWAS_LABEL_NAMES__),
      gwas_label_x_pct=50,
      gwas_label_y_frac=&local_gwas_label_y_frac,
      gwas_label_size=1.8,
      gwas_label_halo_size=1.8,
      gwas_label_angle=0,
      flip1stGWAS_signal=0,
      rm_signals_with_logP_lt=0,
      outputfigname=&batch_output_prefix,
      Use_scaled_pos=1,
      sep_chr_grp=1,
      gwas_sortedby_numchrpos=1,
      angle4xaxis_label=&local_angle4xaxis_label,
      xgrp_y_pos=&local_xgrp_y_pos,
      yoffset_setting=&local_yoffset_setting,
      draw_local_Manhattan=0,
      snp_var=SNP,
      snp_gene_splitter=:,
      target_SNPs=&batch_snp_gene_label,
      LD_SNPs2mark_scatterplot_dots=&local_ld_snps,
      LD_marker_symbol=&local_ld_marker_symbol,
      LD_marker_color=&local_ld_marker_color,
      LD_display_mode=&local_ld_display_mode,
      LD_r2_values=&local_ld_r2_values,
      LD_heatmap_colormodel=&local_ld_heatmap_colors,
      LD_heatmap_legend_title=&local_ld_heatmap_legend_title,
      Keep_order_of_target_SNPs=1
    );
  %end;
%mend;

%render_local_hit_batches;

%macro write_local_hit_html_wrapper;
data _null_;
  file "~/__OUTPUT_PREFIX__.html" lrecl=32767;
  put '<!doctype html>';
  put '<html><head><meta charset="utf-8">';
  put '<title>__HTML_TITLE__</title>';
  put '<style>body{margin:0;padding:16px;font-family:Arial,sans-serif;background:#fff;} img{max-width:100%;height:auto;display:block;}</style>';
  put '</head><body>';
  %do batch_idx=1 %to &n_hit_batches;
    %if &batch_idx=1 %then %do;
      put "<section style=""margin-bottom:24px""><img src=""__OUTPUT_PREFIX__.png"" alt=""__HTML_TITLE__ part &batch_idx of &n_hit_batches""></section>";
    %end;
    %else %do;
      put "<section style=""margin-bottom:24px""><img src=""__OUTPUT_PREFIX___part&batch_idx..png"" alt=""__HTML_TITLE__ part &batch_idx of &n_hit_batches""></section>";
    %end;
  %end;
  put '</body></html>';
run;
%mend;

%write_local_hit_html_wrapper;

ods listing close;
