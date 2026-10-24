options dlcreatedir notes source nodate nonumber;
libname bookout "~/mgxviz_1792866014_2983";
libname bookout clear;
%let BOOK_OUT=%sysget(HOME)/mgxviz_1792866014_2983;
%let P2_OUT=%sysget(HOME)/mgxviz_1792866014_2983;
%let P3_OUT=%sysget(HOME)/mgxviz_1792866014_2983;
%let AF_OUT=%sysget(HOME)/mgxviz_1792866014_2983;
%let AF_REMOTE=mgxviz_1792866014_2983;

/* Zhongshan Cheng. Reusable form of the book's event-glyph oncoprint.
   GPL-3.0. Fixed input schema is documented beside the example.
   Temporary WORK._MGXO_* tables are reserved; caller inputs are never replaced. */
%macro MGXOncoPrint(samples=,genes=,callability=,variants=,
 outlib=work,prefix=mgx_onco,outdir=,image=oncoprint,
 width=7.6,height=5.4,dpi=450,template_file=);
%local i ds did rc field col type bad ns ng maxburden sample_ticks gene_ticks;
%local bad_samples bad_genes bad_call_keys bad_calls bad_events bad_event_classes uncallable_events;
%global MGX_ONCOPRINT_RC;
%let MGX_ONCOPRINT_RC=1;
%if %length(%superq(prefix))=0 or %length(%superq(prefix))>20 %then %do;
 %put ERROR: MGX_ONCOPRINT prefix must be a V7 name of at most 20 characters.;%return;
%end;
%if not %sysfunc(nvalid(%superq(prefix),v7)) or not %sysfunc(nvalid(%superq(image),v7)) %then %do;
 %put ERROR: MGX_ONCOPRINT prefix and image must be V7 names.;%return;
%end;
%if %sysfunc(libref(&outlib)) ne 0 %then %do;
 %put ERROR: MGX_ONCOPRINT outlib is not assigned.;%return;
%end;
%if %length(%superq(outdir))=0 %then %let outdir=%sysfunc(pathname(work));
%if %length(%superq(template_file))=0 %then %do;
 %put ERROR: MGX_ONCOPRINT template_file is required.;%return;
%end;
%let template_file=%sysfunc(dequote(%superq(template_file)));
data _null_;
 length p $2048;p=symget('template_file');
 if substr(p,1,2)='~/' then p=cats(sysget('HOME'),substr(p,2));
 call symputx('template_file',p,'L');
run;
%if not %sysfunc(fileexist(%superq(template_file))) %then %do;
 %put ERROR: MGX_ONCOPRINT cannot read template_file.;%return;
%end;
%if %sysevalf(&width<=0 or &height<=0 or &dpi<=0) %then %do;
 %put ERROR: MGX_ONCOPRINT width height and dpi must be positive.;%return;
%end;
/* Check required column types before any SQL query uses them. */
%do i=1 %to 4;
 %if &i=1 %then %let ds=&samples;
 %else %if &i=2 %then %let ds=&genes;
 %else %if &i=3 %then %let ds=&callability;
 %else %let ds=&variants;
 %if %length(%superq(ds))=0 %then %do;
  %put ERROR: MGX_ONCOPRINT all four input datasets are required.;%return;
 %end;
 %let did=%sysfunc(open(&ds,i));
 %if &did=0 %then %do;
  %put ERROR: MGX_ONCOPRINT cannot open input &ds.;%return;
 %end;
 %if &i=1 %then %let field=sample:C sample_index:N cohort:C;
 %else %if &i=2 %then %let field=gene:C gene_index:N;
 %else %if &i=3 %then %let field=sample:C gene_index:N callable:N;
 %else %let field=sample:C gene:C variant_id:C mutation_class:C;
 %let bad=0;
 %do rc=1 %to %sysfunc(countw(%superq(field),%str( )));
  %let col=%scan(%scan(%superq(field),&rc,%str( )),1,:);
  %let type=%scan(%scan(%superq(field),&rc,%str( )),2,:);
  %let ds=%sysfunc(varnum(&did,&col));
  %if &ds=0 %then %let bad=1;
  %else %if %sysfunc(vartype(&did,&ds)) ne &type %then %let bad=1;
 %end;
 %let rc=%sysfunc(close(&did));
 %if &bad %then %do;
  %put ERROR: MGX_ONCOPRINT required input column missing or wrong type.;%return;
 %end;
%end;
data work._mgxo_samples;set &samples;keep sample sample_index cohort;run;
data work._mgxo_genes;set &genes;keep gene gene_index;run;
data work._mgxo_callability;set &callability;keep sample gene_index callable;run;
data work._mgxo_variants;set &variants;keep sample gene variant_id mutation_class;run;
proc sql noprint;
 select count(*) into :ns trimmed from work._mgxo_samples;
 select count(*) into :ng trimmed from work._mgxo_genes;
 select count(*)-count(distinct sample)+sum(missing(sample) or missing(cohort))
  +count(*)-count(distinct sample_index)
  +sum(sample_index ne int(sample_index) or sample_index<1 or sample_index>&ns)
 into :bad_samples trimmed from work._mgxo_samples;
 select count(*)-count(distinct gene)+sum(missing(gene))
  +count(*)-count(distinct gene_index)
  +sum(gene_index ne int(gene_index) or gene_index<1 or gene_index>&ng)
 into :bad_genes trimmed from work._mgxo_genes;
 create table work._mgxo_bad_keys as select sample,gene_index,count(*) as n
  from work._mgxo_callability group by sample,gene_index having calculated n>1;
 select count(*) into :bad_call_keys trimmed from work._mgxo_bad_keys;
 select count(*) into :bad_calls trimmed from work._mgxo_callability as c
  left join work._mgxo_samples as s on c.sample=s.sample
  left join work._mgxo_genes as g on c.gene_index=g.gene_index
  where missing(s.sample) or missing(g.gene) or c.callable not in (0,1);
 select count(*) into :bad_events trimmed from work._mgxo_variants as v
  left join work._mgxo_samples as s on v.sample=s.sample
  left join work._mgxo_genes as g on v.gene=g.gene
  where missing(s.sample) or missing(g.gene) or missing(v.variant_id)
  or v.mutation_class not in ('Missense','Truncating','Splice','Inframe');
 create table work._mgxo_conflicts as select sample,gene,variant_id,
  count(distinct mutation_class) as n_classes from work._mgxo_variants
  group by sample,gene,variant_id having calculated n_classes>1;
 select count(*) into :bad_event_classes trimmed from work._mgxo_conflicts;
 select count(*) into :uncallable_events trimmed from work._mgxo_variants as v
  left join work._mgxo_genes as g on v.gene=g.gene
  left join work._mgxo_callability as c on v.sample=c.sample and g.gene_index=c.gene_index
  where missing(c.callable) or c.callable ne 1;
quit;
%if %sysevalf(&ns=0 or &ng=0 or &bad_samples>0 or &bad_genes>0 or &bad_call_keys>0 or &bad_calls>0 or &bad_events>0 or &bad_event_classes>0 or &uncallable_events>0,boolean) %then %do;
 %put ERROR: MGX_ONCOPRINT invalid keys indices callability event identities or mutation classes.;
 %return;
%end;
proc sort data=work._mgxo_variants out=work._mgxo_unique nodupkey;
 by sample gene variant_id mutation_class;run;
proc sql;
 create table work._mgxo_counts as select sample,gene,count(*) as n_variants,
  max(mutation_class='Truncating') as truncating,max(mutation_class='Splice') as splice,
  max(mutation_class='Inframe') as inframe,max(mutation_class='Missense') as missense
  from work._mgxo_unique group by sample,gene;
 create table &outlib..&prefix._tiles as
 select s.sample,s.sample_index,s.cohort,g.gene,g.gene_index,c.callable,
  coalesce(m.n_variants,0) as n_variants,coalesce(m.truncating,0) as truncating,
  coalesce(m.splice,0) as splice,coalesce(m.inframe,0) as inframe,coalesce(m.missense,0) as missense
 from work._mgxo_samples as s cross join work._mgxo_genes as g
 left join work._mgxo_callability as c on s.sample=c.sample and g.gene_index=c.gene_index
 left join work._mgxo_counts as m on s.sample=m.sample and g.gene=m.gene;
 create table &outlib..&prefix._frequency as select gene,gene_index,
  sum(callable=1) as n_callable,sum(callable=1 and n_variants>0) as n_mutated,
  case when calculated n_callable>0 then 100*calculated n_mutated/calculated n_callable else . end as pct_mutated
 from &outlib..&prefix._tiles group by gene,gene_index;
 create table work._mgxo_events as select v.*,s.sample_index,g.gene_index,t.n_variants
 from work._mgxo_unique as v inner join work._mgxo_samples as s on v.sample=s.sample
 inner join work._mgxo_genes as g on v.gene=g.gene
 inner join &outlib..&prefix._tiles as t on v.sample=t.sample and v.gene=t.gene;
 create table &outlib..&prefix._burden as select s.sample,s.sample_index as burden_x,
  count(v.variant_id) as n_events,0 as burden_zero from work._mgxo_samples as s
 left join work._mgxo_unique as v on s.sample=v.sample group by s.sample,s.sample_index;
 select max(1,coalesce(max(n_events),0)) into :maxburden trimmed from &outlib..&prefix._burden;
 create table work._mgxo_cohorts as select distinct cohort from work._mgxo_samples order by cohort;
 select sample_index into :sample_ticks separated by ' ' from work._mgxo_samples order by sample_index;
 select gene_index into :gene_ticks separated by ' ' from work._mgxo_genes order by gene_index;
quit;
data &outlib..&prefix._frequency;set &outlib..&prefix._frequency;
 length count_label $32;count_label=cats(n_mutated,'/',n_callable);run;
proc sort data=work._mgxo_events;by sample gene variant_id;run;
data &outlib..&prefix._glyphs;
 set work._mgxo_events;by sample gene;retain event_order;length event_class $24;
 if first.gene then event_order=0;event_order+1;
 event_x=sample_index;event_class=mutation_class;glyph_height=.8/n_variants;
 glyph_low=gene_index-.4+(event_order-1)*glyph_height+.05*glyph_height;
 glyph_high=glyph_low+.9*glyph_height;
run;
data &outlib..&prefix._render;
 set &outlib..&prefix._tiles(in=grid) &outlib..&prefix._frequency(in=bar)
     work._mgxo_samples(in=annotation) &outlib..&prefix._glyphs &outlib..&prefix._burden;
 length background_class $24;
 if grid then do;tile_y=gene_index;
  if callable=1 then background_class='No variant detected';else background_class='Not assessable';end;
 if bar then do;freq_y=gene_index;zero=0;label_x=pct_mutated+3;end;
 if annotation then do;cohort_x=sample_index;strip_y=1;end;
run;
/* CNTLIN preserves arbitrary labels without constructing quoted SAS source. */
data work._mgxo_formats;
 length fmtname $32 type $1 label $256;type='N';
 set work._mgxo_samples(in=s) work._mgxo_genes(in=g);
 if s then do;fmtname='mgxsample';start=sample_index;label=sample;end;
 if g then do;fmtname='mgxgene';start=gene_index;label=gene;end;keep fmtname type start label;
run;
proc format cntlin=work._mgxo_formats;run;
data _null_;
 set work._mgxo_cohorts end=last;
 file "&outdir/&prefix._cohorts.sas" lrecl=32767;
 length color $8 quoted $300;
 array colors[8] $8 _temporary_ ('CX006B8C' 'CXADB5BD' 'CX009E73' 'CXCC79A7' 'CXE69F00' 'CX56B4E9' 'CXD55E00' 'CX332288');
 if _n_=1 then put "discreteattrmap name='cohortcolors';";
 color=colors[mod(_n_-1,8)+1];quoted=quote(trim(cohort),"'");
 if cohort='AML' then color='CX006B8C';
 if cohort='ALL' then color='CXADB5BD';
 put 'value ' quoted ' / fillattrs=(color=' color ');';
 if last then put 'enddiscreteattrmap;';
run;
proc template;%include "&template_file";run;
ods listing gpath="&outdir" image_dpi=&dpi;
ods graphics / reset=index imagename="&image" imagefmt=png width=&width.in height=&height.in;
proc sgrender data=&outlib..&prefix._render template=MGXOncoTemplate;run;
ods listing gpath="&outdir" image_dpi=96;
ods graphics / reset=index imagename="&image" imagefmt=svg width=&width.in height=&height.in;
proc sgrender data=&outlib..&prefix._render template=MGXOncoTemplate;run;
%let MGX_ONCOPRINT_RC=0;
%put NOTE: MGX_ONCOPRINT_COMPLETE samples=&ns genes=&ng.;
%mend;


ods listing gpath="&BOOK_OUT";

/* LESSON_BEGIN */
/* Invented 2x2 assay: duplicate event, multiple classes, 100% and 0/0. */
data work.samples;length sample $3 cohort $8;
 input sample $ sample_index cohort $;datalines;
C1 1 ArmA
C2 2 ArmB
;
run;
data work.genes;length gene $8;input gene $ gene_index;datalines;
GENE1 1
GENE2 2
;
run;
data work.calls;length sample $3;input sample $ gene_index callable;datalines;
C1 1 1
C2 1 1
C1 2 0
;
run;
data work.variants;length sample $3 gene $8 variant_id $10 mutation_class $12;
 input sample $ gene $ variant_id $ mutation_class $;datalines;
C1 GENE1 e1 Missense
C1 GENE1 e1 Missense
C1 GENE1 e2 Truncating
C2 GENE1 e3 Splice
;
run;
%MGXOncoPrint(samples=work.samples,genes=work.genes,
 callability=work.calls,variants=work.variants,
 prefix=edge,outdir=&BOOK_OUT,image=oncoprint_edges,width=6,height=4,
 template_file="~/mgxviz_1792866014_2983_oncoprint_template.sas");
data work.edge_status;status=&MGX_ONCOPRINT_RC;run;

/* LESSON_END */

/* Reader capture: export values without changing the analysis. */
data work.reader_export;set work.edge_tiles;format _numeric_ best32.;run;
proc export data=work.reader_export outfile="~/mgxviz_1792866014_2983/edge_tiles.csv" dbms=csv replace;run;
data work.reader_export;set work.edge_frequency;format _numeric_ best32.;run;
proc export data=work.reader_export outfile="~/mgxviz_1792866014_2983/edge_frequency.csv" dbms=csv replace;run;
data work.reader_export;set work.edge_burden;format _numeric_ best32.;run;
proc export data=work.reader_export outfile="~/mgxviz_1792866014_2983/edge_burden.csv" dbms=csv replace;run;
data work.reader_export;set work.edge_glyphs;format _numeric_ best32.;run;
proc export data=work.reader_export outfile="~/mgxviz_1792866014_2983/edge_glyphs.csv" dbms=csv replace;run;
data work.reader_export;set work.edge_render;format _numeric_ best32.;run;
proc export data=work.reader_export outfile="~/mgxviz_1792866014_2983/edge_render.csv" dbms=csv replace;run;
data work.reader_export;set work.edge_status;format _numeric_ best32.;run;
proc export data=work.reader_export outfile="~/mgxviz_1792866014_2983/edge_status.csv" dbms=csv replace;run;
data work.reader_environment;length sas_version $64 platform $32 process_id $40; sas_version="&sysvlong";platform="&sysscp";process_id="&sysprocessid";run;
proc export data=work.reader_environment outfile="~/mgxviz_1792866014_2983/environment.csv" dbms=csv replace;run;
