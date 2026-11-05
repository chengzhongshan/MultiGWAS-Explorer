/* Parameters real_input, real_expected_n and real_metadata_dir are set by
   run_real_gtex_boxplots.sh. Input is a small, real two-gene matrix subset.
   No full-matrix download is performed. The metadata library is read-only. */
libname gtexmeta "&real_metadata_dir" access=readonly;
filename realraw "&real_input";
data work.real_expression;
  length rownames $64 sample_id $200;
  infile realraw dlm='09'x dsd firstobs=2 truncover lrecl=4096;
  input rownames :$64. sample_id :$200. raw_exp :best32.;
run;
/* Delete only the uniquely named upload created for this invocation. */
data _null_;
  rc=fdelete('realraw');
  if rc ne 0 then put 'WARNING: Could not remove this run''s temporary GTEx upload.';
run;
filename realraw clear;

%macro plot_real_gtex;
  %local input_n sample_n metadata_n metadata_unique bad_n joined_n plotted_n genes_n;
  %if not %sysfunc(exist(gtexmeta.umap)) %then %do;
    %put ERROR: GTEX_REAL: The metadata library must contain UMAP with seq_id cluster sex Tissue.;
    %return;
  %end;
  proc sql noprint;
    select count(*),count(distinct sample_id) into :input_n trimmed,:sample_n trimmed
      from work.real_expression;
    select count(*),count(distinct seq_id) into :metadata_n trimmed,:metadata_unique trimmed
      from gtexmeta.umap;
  quit;
  %if &input_n ne &real_expected_n or &metadata_n ne &metadata_unique %then %do;
    %put ERROR: GTEX_REAL: Input row count differs from the local conversion or phenotype sample IDs are duplicated.;
    %return;
  %end;
  proc sql;
    create table work.real_joined as
    select a.*,b.seq_id,b.cluster,b.sex,b.Tissue
      from work.real_expression as a left join gtexmeta.umap as b
      on a.sample_id=b.seq_id;
  quit;
  proc sql noprint;
    select count(*) into :joined_n trimmed from work.real_joined;
    select count(*) into :bad_n trimmed from work.real_joined
      where missing(seq_id) or missing(cluster) or sex not in ('F','M')
        or missing(raw_exp) or raw_exp<0;
  quit;
  %if &joined_n ne &input_n or &bad_n ne 0 %then %do;
    %put ERROR: GTEX_REAL: Incomplete or ambiguous sample/phenotype join or invalid expression values.;
    %return;
  %end;
  data work.real_gtex_cache;
    set work.real_joined;
    /* Match the original boxplot macro's transformation and outlier rule. */
    exp=log2(raw_exp+1);
    label exp='log2(normalized expression + 1)';
    if index(upcase(cluster),'OUTLIER')=0;
  run;
  %Boxplots4GenesInGTExV8_Cached(
    genes=MAP3K19 R3HDM1,
    dsdout=work.real_gtex_plot,
    bysex=1,
    UseGeneratedDsd=1,
    PreviousDsd=real_gtex_cache,
    Lib4PreviousDsd=work,
    WhereFilters4Boxplot=%str(),
    boxplot_width=2800,
    boxplot_height=1600
  );
  %if not %sysfunc(exist(work.real_gtex_plot)) %then %return;
  proc sql noprint;
    select count(*),count(distinct rownames) into :plotted_n trimmed,:genes_n trimmed
      from work.real_gtex_plot;
  quit;
  %if &plotted_n=0 or &genes_n ne 2 %then %do;
    %put ERROR: GTEX_REAL: The plot must contain both requested genes.;
    %return;
  %end;
  title 'Real GTEx V8 expression: sample and tissue counts';
  proc sql;
    select rownames,count(*) as observations,count(distinct seq_id) as samples,
      count(distinct cluster) as tissues from work.real_gtex_plot group by rownames;
  quit;
  title 'Real GTEx V8: log2(normalized expression + 1) by tissue and sex';
  proc means data=work.real_gtex_plot n mean median min max;
    class rownames cluster sex;
    var exp;
  run;
  title;
  %put GTEX_REAL_PASS: input=&input_n samples=&sample_n joined=&joined_n plotted=&plotted_n genes=&genes_n;
%mend;
%plot_real_gtex;
libname gtexmeta clear;
