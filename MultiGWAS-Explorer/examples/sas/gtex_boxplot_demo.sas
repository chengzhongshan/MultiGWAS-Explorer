/* Synthetic learning example: these are not measured GTEx expression values.
   Run after loading Boxplots4GenesInGTExV8_Cached.sas in the same submission. */
data work.gtex_demo_cache;
  length rownames $16 cluster $12 sex $1;
  do gene=1 to 3;
    if gene=1 then rownames='MAP3K19';
    else if gene=2 then rownames='R3HDM1';
    else rownames='OTHER';
    do tissue=1 to 3;
      if tissue=1 then cluster='Liver';
      else if tissue=2 then cluster='Lung';
      else cluster='Whole_Blood';
      do sample=1 to 8;
        sex='F'; exp=gene*2+tissue+sample/4; output;
        sex='M'; exp=gene*2+tissue+sample/4+0.5; output;
      end;
    end;
  end;
  keep rownames cluster sex sample exp;
run;

%Boxplots4GenesInGTExV8_Cached(
  genes=MAP3K19 R3HDM1,
  dsdout=work.gtex_demo_selected,
  bysex=1,
  UseGeneratedDsd=1,
  PreviousDsd=gtex_demo_cache,
  Lib4PreviousDsd=work,
  WhereFilters4Boxplot=%str(),
  boxplot_width=1200,
  boxplot_height=1000
);

/* Print a few selected observations and the numbers behind the boxplots. */
title 'Synthetic example: first 12 selected observations';
proc print data=work.gtex_demo_selected(obs=12) noobs;
  var rownames cluster sex sample exp;
run;
title 'Synthetic example: expression summary by gene, tissue and sex';
proc means data=work.gtex_demo_selected n mean median min max;
  class rownames cluster sex;
  var exp;
run;
title;

proc sql noprint;
  select count(*) into :demo_source_n trimmed from work.gtex_demo_cache;
  select count(*), count(distinct rownames) into :demo_selected_n trimmed,
    :demo_gene_n trimmed from work.gtex_demo_selected;
  select count(*) into :demo_other_n trimmed from work.gtex_demo_selected
    where rownames not in ('MAP3K19','R3HDM1');
quit;
data _null_;
  if &demo_source_n=144 and &demo_selected_n=96 and
     &demo_gene_n=2 and &demo_other_n=0 then
    put 'GTEX_DEMO_PASS: source=144 selected=96 genes=2';
  else put 'ERROR: GTEx demo subset validation failed.';
run;
