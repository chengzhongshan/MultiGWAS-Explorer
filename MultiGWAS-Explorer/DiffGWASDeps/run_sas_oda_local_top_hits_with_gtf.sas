/*
Run this in SAS ODA after uploading:
  1) a wide differential GWAS subset .tsv.gz with the expected beta / SE / P columns

This script follows the user's GetTopHits_and_Make_Local_Manhattan.sas workflow
more closely by:
  - selecting top loci from the wide differential GWAS subset
  - annotating those loci with HaploReg genes
  - building local windows around each hit
  - computing Z-score variables from beta / SE
  - calling SNP_Local_Manhattan_With_GTF to draw local Manhattan plots with
    gene tracks

Note:
  SNP_Local_Manhattan_With_GTF is expected to be available automatically inside
  the SAS ODA environment, consistent with the user's existing workflow.
*/

*options mprint mlogic symbolgen;

%let top_hit_focus_pvar=__TOP_HIT_FOCUS_PVAR__;
%let top_hit_mode=__TOP_HIT_MODE__;
%let top_hit_filter_expr=__TOP_HIT_FILTER_EXPR__;
%let top_hit_signal_thrshd=__TOP_HIT_SIGNAL_THRSHD__;
%let top_hit_signal_thrshds=__TOP_HIT_SIGNAL_THRSHDS__;
%let top_hit_dist_bp=__TOP_HIT_DIST_BP__;
%let local_window_bp=__LOCAL_WINDOW_BP__;
%let local_max_hits_per_fig=__LOCAL_MAX_HITS_PER_FIG__;
%let local_top_hits_csv_basename=__LOCAL_TOP_HITS_CSV_BASENAME__;
%let lth_input_csv=__LOCAL_TOP_HITS_INPUT_CSV_BASENAME__;
%let target_snp_list=__TARGET_SNP_LIST__;
%let target_snp_gene_map=__TARGET_SNP_GENES__;
%let common_assoc_pvars=__COMMON_ASSOC_P_VARS__;
%let prep_only=__PREP_ONLY__;

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
%let gtf_yaxis_offset4max=__GTF_YAXIS_OFFSET4MAX__;
%let gtf_yoffset4textlabels=__GTF_YOFFSET4TEXTLABELS__;
%let gtf_yoffset4maxdrawmarkersontop=__GTF_YOFFSET4MAX_DRAWMARKERSONTOP__;
%let gtf_label_snps=__GTF_LABEL_SNPS__;
%let gtf_label_text_rotate_angle=__GTF_LABEL_TEXT_ROTATE_ANGLE__;
%let gtf_include_non_protein_coding=__GTF_INCLUDE_NON_PROTEIN_CODING__;

/* Start each local-GTF batch from a clean WORK library to reduce
   cross-batch accumulation in persistent SAS ODA sessions. */
proc datasets library=work kill nolist memtype=(data view catalog);
quit;

ods _all_ close;
ods html5 file="~/__OUTPUT_HTML__"
  options(bitmap_mode='inline')
  style=HTMLBlue;
ods graphics on / outputfmt=png;

%if %length(&fm_libpath) > 0 %then %do;
  libname FM "&fm_libpath";
%end;

%include "~/get_top_signal_within_dist.sas";
%include "~/__GET_GTF_MACRO_BASENAME__";
%include "~/adj_grpnum4close_gene_bed_regs.sas";
%include "~/Multgscatter_with_gene_exons.sas";
%include "~/map_grp_assoc2gene4covidsexgwas.sas";
%include "~/SNP_Local_Manhattan_With_GTF.sas";
%include "~/Lattice_gscatter_over_bed_track.sas";

__GTF_IMPORT_BLOCK__

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

  %if %sysevalf(%superq(lth_input_csv)=,boolean) %then %return;

  filename _reqcsv "~/&lth_input_csv";
  %if not %sysfunc(fexist(_reqcsv)) %then %do;
    %put NOTE: Requested top-hit CSV ~/&lth_input_csv was not found in SAS ODA home. Using computed top hits instead.;
    filename _reqcsv clear;
    %return;
  %end;

  proc import datafile="~/&lth_input_csv"
    out=&outdsd
    dbms=csv
    replace;
    guessingrows=max;
  run;
  filename _reqcsv clear;

  %if not %sysfunc(exist(&outdsd)) %then %return;

  %_find_first_column(lib=WORK,mem=%upcase(&outdsd),outvar=req_chr_var,candidates=CHR chromosome chrom);
  %_find_first_column(lib=WORK,mem=%upcase(&outdsd),outvar=req_bp_var,candidates=BP POS position);
  %_find_first_column(lib=WORK,mem=%upcase(&outdsd),outvar=req_snp_var,candidates=SNP rsid markername MarkerName);
  %_find_first_column(lib=WORK,mem=%upcase(&outdsd),outvar=req_hit_order_var,candidates=hit_order panel_order order);

  %if %superq(req_chr_var)= or %superq(req_bp_var)= or %superq(req_snp_var)= %then %do;
    %put WARNING: Requested top-hit CSV lacks resolvable CHR/BP/SNP columns, so it will be ignored.;
    %put WARNING: Resolved columns: chr=&req_chr_var bp=&req_bp_var snp=&req_snp_var hit_order=&req_hit_order_var;
    %return;
  %end;

  data &outdsd;
    length SNP $128 _chr_text _bp_text _snp_text _hit_order_text $256;
    set &outdsd;
    _chr_text=strip(vvaluex("&req_chr_var"));
    _bp_text=strip(vvaluex("&req_bp_var"));
    _snp_text=strip(vvaluex("&req_snp_var"));
    %if %superq(req_hit_order_var) ne %then %do;
    _hit_order_text=strip(vvaluex("&req_hit_order_var"));
    %end;
    %else %do;
    _hit_order_text='';
    %end;

    if upcase(_chr_text)='X' then CHR=23;
    else if upcase(_chr_text)='Y' then CHR=24;
    else if upcase(_chr_text) in ('M','MT') then CHR=25;
    else CHR=input(_chr_text,best32.);
    BP=input(_bp_text,best32.);
    SNP=_snp_text;
    hit_order=input(_hit_order_text,best32.);
    if missing(hit_order) then hit_order=_n_;
    if missing(CHR) or missing(BP) or missing(SNP) then delete;
    keep hit_order CHR BP SNP;
  run;

  proc sort data=&outdsd nodupkey;
    by hit_order CHR BP SNP;
  run;

  %let requested_top_hits_loaded=1;
  %put NOTE: Loaded requested top-hit CSV ~/&lth_input_csv for targeted local GTF plotting.;
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
  %global effective_gtf_dsd effective_gtf_dist2snp effective_gtf_pct4neg_y effective_gtf_dist2sep_genes effective_gtf_grp_font_size;
  %local _source_gtf_dsd;
  %local _n_genes_for_plot _n_exons_for_plot _n_non_protein_for_plot _n_tx_for_plot;
  %let effective_gtf_dist2snp=&gtf_dist2snp;
  %let effective_gtf_pct4neg_y=&gtf_pct4neg_y;
  %let effective_gtf_dist2sep_genes=&gtf_dist2sep_genes;
  %let effective_gtf_grp_font_size=10;
  %if %sysevalf(&local_window_bp > &effective_gtf_dist2snp) %then %do;
    %put NOTE: Expanding local GTF half-window from &effective_gtf_dist2snp bp to the requested SNP-centered local window (&local_window_bp bp).;
    %let effective_gtf_dist2snp=&local_window_bp;
  %end;

  %if %sysfunc(exist(&gtf_local_dsd)) %then %do;
    %let _source_gtf_dsd=&gtf_local_dsd;
    %let effective_gtf_dsd=&gtf_local_dsd;
    %put NOTE: Reusing uploaded local GTF subset &gtf_local_dsd for local gene-track plotting.;
  %end;
  %else %do;
    %let effective_gtf_dsd=&gtf_local_dsd._for_plot;
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
    %let _source_gtf_dsd=&gtf_local_dsd;
 
    data &gtf_local_dsd._for_plot;
      length chr 8 chr_text $64 ensembl $64 type $32 genesymbol gene $256 protein_coding 8 original_protein_coding 8 _bio_text $32767;
      set &_source_gtf_dsd(rename=(chr=chr_text));
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
  %end;

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

  %put NOTE: Using locally generated rich GTF dataset &effective_gtf_dsd for local gene-track plotting.;
  %put NOTE: Local GTF stats: genes=&_n_genes_for_plot transcripts=&_n_tx_for_plot exons=&_n_exons_for_plot non_protein_coding_features=&_n_non_protein_for_plot include_non_protein_coding=&gtf_include_non_protein_coding display_window_bp=&effective_gtf_dist2snp.;
%mend;

%macro _auto_tune_gene_track_ratio(top_hits_dsd=,gtf_dsd=);
  %global effective_gtf_pct4neg_y effective_gtf_dist2sep_genes effective_gtf_grp_font_size;
  %local _max_gene_count_per_locus _recommended_pct4neg_y _recommended_dist2sep_genes _recommended_grp_font_size _display_window_bp _normalized_base_dist2sep_genes;
  %let _max_gene_count_per_locus=0;
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
    create table _gtf_gene_density as
    select a.SNP,
           count(distinct coalescec(b.genesymbol,b.gene)) as gene_count
    from &top_hits_dsd as a
    left join &gtf_dsd as b
      on a.CHR=b.chr
     and b.type='gene'
     and b.st <= (a.BP + &effective_gtf_dist2snp)
     and b.en >= (a.BP - &effective_gtf_dist2snp)
    group by a.SNP
    ;

    select max(gene_count)
      into :_max_gene_count_per_locus trimmed
    from _gtf_gene_density
    ;
  quit;

  %if %sysevalf(%superq(_max_gene_count_per_locus)=,boolean) %then %let _max_gene_count_per_locus=0;

  %if %sysevalf(&_max_gene_count_per_locus > 120) %then %do;
    %let _recommended_pct4neg_y=%sysfunc(max(&gtf_pct4neg_y,4.20));
    %let _recommended_dist2sep_genes=%sysfunc(max(&_normalized_base_dist2sep_genes,%sysfunc(max(250000,%sysfunc(int(%sysevalf(&_display_window_bp*0.020)))))));
    %let _recommended_grp_font_size=6;
  %end;
  %else %if %sysevalf(&_max_gene_count_per_locus > 80) %then %do;
    %let _recommended_pct4neg_y=%sysfunc(max(&gtf_pct4neg_y,3.40));
    %let _recommended_dist2sep_genes=%sysfunc(max(&_normalized_base_dist2sep_genes,%sysfunc(max(200000,%sysfunc(int(%sysevalf(&_display_window_bp*0.015)))))));
    %let _recommended_grp_font_size=7;
  %end;
  %else %if %sysevalf(&_max_gene_count_per_locus > 50) %then %do;
    %let _recommended_pct4neg_y=%sysfunc(max(&gtf_pct4neg_y,2.60));
    %let _recommended_dist2sep_genes=%sysfunc(max(&_normalized_base_dist2sep_genes,%sysfunc(max(150000,%sysfunc(int(%sysevalf(&_display_window_bp*0.010)))))));
    %let _recommended_grp_font_size=7;
  %end;
  %else %if %sysevalf(&_max_gene_count_per_locus > 30) %then %do;
    %let _recommended_pct4neg_y=%sysfunc(max(&gtf_pct4neg_y,2.00));
    %let _recommended_dist2sep_genes=%sysfunc(max(&_normalized_base_dist2sep_genes,%sysfunc(max(125000,%sysfunc(int(%sysevalf(&_display_window_bp*0.008)))))));
    %let _recommended_grp_font_size=8;
  %end;
  %else %if %sysevalf(&_max_gene_count_per_locus > 15) %then %do;
    %let _recommended_pct4neg_y=%sysfunc(max(&gtf_pct4neg_y,1.50));
    %let _recommended_dist2sep_genes=%sysfunc(max(&_normalized_base_dist2sep_genes,%sysfunc(max(100000,%sysfunc(int(%sysevalf(&_display_window_bp*0.006)))))));
    %let _recommended_grp_font_size=9;
  %end;
  %else %if %sysevalf(&_max_gene_count_per_locus > 8) %then %do;
    %let _recommended_pct4neg_y=%sysfunc(max(&gtf_pct4neg_y,1.20));
    %let _recommended_dist2sep_genes=%sysfunc(max(&_normalized_base_dist2sep_genes,%sysfunc(max(80000,%sysfunc(int(%sysevalf(&_display_window_bp*0.004)))))));
    %let _recommended_grp_font_size=10;
  %end;
  %else %if %sysevalf(&_max_gene_count_per_locus > 0) %then %do;
    %let _recommended_pct4neg_y=%sysfunc(max(&gtf_pct4neg_y,1.10));
    %let _recommended_dist2sep_genes=%sysfunc(max(&_normalized_base_dist2sep_genes,%sysfunc(max(60000,%sysfunc(int(%sysevalf(&_display_window_bp*0.003)))))));
    %let _recommended_grp_font_size=10;
  %end;

  %let effective_gtf_pct4neg_y=&_recommended_pct4neg_y;
  %let effective_gtf_dist2sep_genes=&_recommended_dist2sep_genes;
  %let effective_gtf_grp_font_size=&_recommended_grp_font_size;
  %put NOTE: Auto-tuned local GTF layout using max_genes_per_locus=&_max_gene_count_per_locus pct4neg_y=&effective_gtf_pct4neg_y (base=&gtf_pct4neg_y) dist2sep_genes=&effective_gtf_dist2sep_genes (base=&gtf_dist2sep_genes) grp_font_size=&effective_gtf_grp_font_size.;
%mend;

%macro _set_force_signal_xaxis_bounds(signal_dsd=,top_hits_dsd=);
  %global force_lattice_xaxis_viewmin force_lattice_xaxis_viewmax;
  %local _n_target_snps _n_target_chrs _center_bp _target_min_bp _target_max_bp;
  %let force_lattice_xaxis_viewmin=;
  %let force_lattice_xaxis_viewmax=;

  proc sql noprint;
    select count(distinct strip(SNP)),
           count(distinct CHR)
      into :_n_target_snps trimmed,
           :_n_target_chrs trimmed
    from &top_hits_dsd
    where prxmatch('/^rs[0-9]+$/i', strip(SNP)) > 0
    ;
  quit;

  %if %sysevalf(%superq(_n_target_snps)=,boolean) %then %let _n_target_snps=0;
  %if %sysevalf(%superq(_n_target_chrs)=,boolean) %then %let _n_target_chrs=0;

  %if %eval(&_n_target_snps=1 and &_n_target_chrs=1) %then %do;
    proc sql noprint;
      select int(mean(BP))
        into :_center_bp trimmed
      from &top_hits_dsd
      where not missing(BP)
      ;
    quit;

    %if %sysevalf(%superq(_center_bp)^=,boolean) %then %do;
      %let _target_min_bp=%sysfunc(int(%sysevalf(&_center_bp-&effective_gtf_dist2snp)));
      %let _target_max_bp=%sysfunc(int(%sysevalf(&_center_bp+&effective_gtf_dist2snp)));
      %if %sysevalf(&_target_min_bp<1) %then %let _target_min_bp=1;
      %let force_lattice_xaxis_viewmin=&_target_min_bp;
      %let force_lattice_xaxis_viewmax=&_target_max_bp;
      %put NOTE: Forcing the displayed x-axis to the requested SNP-centered window [&force_lattice_xaxis_viewmin, &force_lattice_xaxis_viewmax].;
    %end;
  %end;
  %else %do;
    %put NOTE: Skipping forced x-axis bounds because the local GTF plot includes &_n_target_snps target SNPs across &_n_target_chrs chromosomes.;
  %end;
%mend;

%macro _maybe_expand_gtf_dist2snp(top_hits_dsd=,gtf_dsd=);
  %local _gtf_lib _gtf_mem _nearby_gene_count;
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

  %if %superq(gtf_chr_var)= or %superq(gtf_start_var)= or %superq(gtf_end_var)= or %superq(gtf_gene_var)= %then %return;

  %let _nearby_gene_count=0;
  proc sql noprint;
    select count(*) into: _nearby_gene_count trimmed
    from &top_hits_dsd as a, &gtf_dsd as b
    where not missing(strip(cats(b.&gtf_gene_var)))
      and b.&gtf_start_var <= b.&gtf_end_var
      and upcase(compress(strip(cats(a.CHR)),'CHR')) = upcase(compress(strip(cats(b.&gtf_chr_var)),'CHR'))
      and (a.BP between (b.&gtf_start_var-&gtf_dist2snp) and (b.&gtf_end_var+&gtf_dist2snp))
      %if %superq(gtf_feature_var) ne %then %do;
      and (missing(strip(cats(b.&gtf_feature_var))) or lowcase(strip(cats(b.&gtf_feature_var)))='gene')
      %end;
    ;
  quit;

  %if %sysevalf(&_nearby_gene_count=0) %then %do;
    %let effective_gtf_dist2snp=&local_window_bp;
    %put NOTE: No genes were found within the default GTF distance of &gtf_dist2snp bp. Expanding gene-track search to the local window size (&effective_gtf_dist2snp bp).;
  %end;
%mend;

__WIDE_IMPORT_BLOCK__

proc sort data=scz_mh;
  by CHR BP;
run;

%if %upcase(&top_hit_mode)=COMMON_ASSOCIATION %then %do;
data scz_mh;
  set scz_mh;
  COMMON_ASSOC_P=min(of &common_assoc_pvars);
run;
%end;

%_load_requested_target_snps(outdsd=requested_target_snps);

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
%end;
%else %do;
data top_hit_candidates;
  set scz_mh;
  if &top_hit_filter_expr;
run;

%macro _pick_top_hits_by_thr;
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
      from top_hit_candidates
      where (&top_hit_focus_pvar>0) and (&top_hit_focus_pvar<&_thr);
    quit;
    %if %sysevalf(&_n_hits_this>0) %then %do;
      %let top_hit_signal_thrshd=&_thr;
      %goto _picked;
    %end;
  %end;
  %_picked:
%mend;
%_pick_top_hits_by_thr;

%get_top_signal_within_dist(
  dsdin=top_hit_candidates,
  grp_var=CHR,
  signal_var=&top_hit_focus_pvar,
  select_smallest_signal=1,
  pos_var=BP,
  pos_dist_thrshd=&top_hit_dist_bp,
  dsdout=top_hit4diffp_raw,
  signal_thrshd=&top_hit_signal_thrshd
);
%end;

%if not %sysevalf(&requested_target_snps_loaded,boolean) %then %do;
data top_hit4diffp_raw;
  set top_hit4diffp_raw;
  requested_hit_order=.;
run;
%end;

%_load_requested_target_snp_genes(outdsd=requested_target_snp_genes);
%if not %sysfunc(exist(work.requested_target_snp_genes)) %then %do;
  data requested_target_snp_genes;
    length SNP $128 gene $256;
    stop;
  run;
%end;

%_load_req_top_hits_csv(outdsd=requested_top_hits_csv);

%if %sysevalf(&requested_top_hits_loaded,boolean) %then %do;
proc sql;
  create table top_hit4diffp_raw as
  select a.*,
         b.hit_order as requested_hit_order
  from scz_mh as a
  inner join requested_top_hits_csv as b
    on a.CHR=b.CHR
   and a.BP=b.BP
   and a.SNP=b.SNP
  order by b.hit_order, a.CHR, a.BP, a.SNP
  ;
quit;
%put NOTE: Requested local-top-hit CSV is being treated as the explicit source of truth for this GTF batch, so previously selected loci will not be re-pruned by distance or threshold filtering.;
%end;

proc sql noprint;
  select SNP into: top_snps separated by ' '
  from top_hit4diffp_raw
  order by coalesce(requested_hit_order, 999999999), CHR, BP;
quit;

%if %superq(top_snps)= %then %do;
  %put ERROR: No top SNPs passed the local top-hit filter and threshold settings.;
  %abort 255;
%end;

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

%_ensure_effective_gtf_dsd(top_hits_dsd=top_hit4diffp_raw);

%_prepare_gtf_gene_fallback(
  top_hits_dsd=top_hit4diffp_raw,
  gtf_dsd=&effective_gtf_dsd,
  outdsd=snps2genes_gtf_fallback
);

proc sql;
  create table top_hit4diffp as
  select a.*,
         coalescec(u.gene, b.gene, c.gtf_gene) as gene length=256,
         case
           when not missing(u.gene) then 'USER'
           when not missing(b.gene) then 'HaploReg'
           when not missing(c.gtf_gene) then 'GTF'
           else 'NA'
         end as gene_source length=16,
         catx(':', a.SNP, coalescec(u.gene, b.gene, c.gtf_gene, 'NA')) as snp_gene length=128
  from top_hit4diffp_raw as a
  left join requested_target_snp_genes as u
    on upcase(strip(a.SNP))=upcase(strip(u.SNP))
  left join snps2genes_clean as b
    on a.SNP=b.rsid
  left join snps2genes_gtf_fallback as c
    on a.SNP=c.rsid
  ;

  create table top_hit_groups as
  select distinct
         t.snp_gene,
         coalesce(
         %if %sysevalf(&requested_top_hits_loaded,boolean) %then %do;
         r.hit_order,
         %end;
         t.requested_hit_order,
         .
         )
           as hit_order
  from top_hit4diffp as t
  %if %sysevalf(&requested_top_hits_loaded,boolean) %then %do;
  left join requested_top_hits_csv as r
    on t.CHR=r.CHR
   and t.BP=r.BP
   and t.SNP=r.SNP
  %end;
  order by calculated hit_order, t.CHR, t.BP, t.SNP
  ;

  create table top_hit4diffp_gtf_ready as
  select *
  from top_hit4diffp
      where prxmatch('/^rs[0-9]+$/i', strip(SNP)) > 0
      ;

  create table top_local_signals as
  select a.*, catx(':', b.SNP, coalescec(b.gene,'NA')) as snp_gene length=128
  from scz_mh as a, top_hit4diffp_gtf_ready as b
  where a.CHR=b.CHR
    and a.BP between (b.BP-&local_window_bp) and (b.BP+&local_window_bp)
  ;
quit;

data top_hit_groups;
  set top_hit_groups;
  hit_order=_n_;
  panel_index=ceil(hit_order / max(1,&local_max_hits_per_fig));
run;

%_maybe_expand_gtf_dist2snp(
  top_hits_dsd=top_hit4diffp_gtf_ready,
  gtf_dsd=&effective_gtf_dsd
);
%_auto_tune_gene_track_ratio(
  top_hits_dsd=top_hit4diffp_gtf_ready,
  gtf_dsd=&effective_gtf_dsd
);
%_set_force_signal_xaxis_bounds(
  signal_dsd=top_local_signals,
  top_hits_dsd=top_hit4diffp_gtf_ready
);

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
         g.panel_index,
         t.CHR,
         t.BP,
         t.SNP,
         t.A1 as EFFECT_ALLELE length=8,
         t.A2 as OTHER_ALLELE length=8,
         t.A2 as REFERENCE_ALLELE length=8,
         t.A1 as ALTERNATIVE_ALLELE length=8,
         t.gene,
         t.snp_gene,
         t.&top_hit_focus_pvar as focus_signal,
         %if %superq(top_hit_export_extra_vars) ne %then %do;
         &top_hit_export_extra_vars
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

%macro _emit_or_plot_local_gtf;
  %global top_snps effective_gtf_label_snps;
  %if %sysevalf(&prep_only,boolean) %then %do;
    data _null_;
      file "~/__OUTPUT_HTML__" lrecl=32767;
      put '<!doctype html>';
      put '<html><head><meta charset="utf-8"><title>Local Top-Hit GTF Prep</title></head>';
      put '<body style="font-family:Arial,Helvetica,sans-serif;margin:24px">';
      put '<h1 style="font-size:20px">Local top-hit GTF prep completed</h1>';
      put "<p>Exported top-hit CSV: &local_top_hits_csv_basename</p>";
      put '</body></html>';
    run;
    ods html5 close;
    %put NOTE: Prep-only mode complete; exported local top-hit CSV and skipped GTF plotting.;
    %return;
  %end;

  proc sql noprint;
    select SNP into: top_snps separated by ' '
    from top_hit4diffp_gtf_ready as t
    inner join top_hit_groups as g
      on t.snp_gene=g.snp_gene
    order by g.hit_order, t.CHR, t.BP;
  quit;

  %if %superq(top_snps)= %then %do;
    %put ERROR: No rsID-style top SNPs remain after filtering invalid local GTF targets.;
    data _null_;
      file "~/__OUTPUT_HTML__" lrecl=32767;
      put '<!doctype html>';
      put '<html><head><meta charset="utf-8"><title>Local Top-Hit GTF Error</title></head>';
      put '<body style="font-family:Arial,Helvetica,sans-serif;margin:24px">';
      put '<h1 style="font-size:20px">Local top-hit GTF plotting failed</h1>';
      put '<p>No rsID-style top SNPs remained after filtering invalid local GTF targets.</p>';
      put '</body></html>';
    run;
    ods html5 close;
    %return;
  %end;

  %let effective_gtf_label_snps=&top_snps;
  %if %superq(gtf_label_snps) ne %then %let effective_gtf_label_snps=&gtf_label_snps;

  %SNP_Local_Manhattan_With_GTF(
    gwas_dsd=top_local_signals,
    chr_var=CHR,
    AssocPVars=&gtf_assoc_pvars,
    SNP_IDs=&top_snps,
    SNP_Var=SNP,
    Pos_Var=BP,
    gtf_dsd=&effective_gtf_dsd,
    dist2snp=&effective_gtf_dist2snp,
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
    SNPs2label_scatterplot_dots=&effective_gtf_label_snps,
    %if %superq(gtf_label_text_rotate_angle) ne %then %do;
    text_rotate_angle=&gtf_label_text_rotate_angle,
    %end;
    %if %superq(gtf_yaxis_offset4max) ne %then %do;
    yaxis_offset4max=&gtf_yaxis_offset4max,
    %end;
    yoffset4max_drawmarkersontop=&gtf_yoffset4maxdrawmarkersontop,
    Yoffset4textlabels=&gtf_yoffset4textlabels,
    verbose=0
  );

  ods html5 close;
%mend;
%_emit_or_plot_local_gtf;
