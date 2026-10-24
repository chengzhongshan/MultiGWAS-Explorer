data work.samples;
 length sample $3 cohort $8;
 input sample $ cohort $;
 sample_index=input(substr(sample,2),best.);
 datalines;
S01 AML
S02 AML
S03 AML
S04 AML
S05 AML
S06 AML
S07 ALL
S08 ALL
S09 ALL
S10 ALL
S11 ALL
S12 ALL
;
run;

data work.genes;
 length gene $8;
 input gene $ gene_index;
 datalines;
TP53 8
DNMT3A 7
RUNX1 6
FLT3 5
NPM1 4
TET2 3
IDH1 2
ASXL1 1
;
run;

/* CALLABLE=1 means assessable; zero means unavailable in this fixture. */
data work.callability;
 set work.samples;
 do gene_index=1 to 8;
   callable=1;
   if (sample='S04' and gene_index=8) or
      (sample='S10' and gene_index=2) then callable=0;
   output;
 end;
 keep sample gene_index callable;
run;

/* Event table preserves individual records before collapsing. */
data work.variants;
 length sample $3 gene $8 variant_id $10 mutation_class $12;
 input sample $ gene $ variant_id $ mutation_class $;
 datalines;
S01 TP53 v001 Truncating
S02 TP53 v002 Missense
S03 TP53 v003 Missense
S03 TP53 v004 Truncating
S05 TP53 v005 Splice
S07 TP53 v006 Missense
S08 TP53 v007 Truncating
S11 TP53 v008 Missense
S01 DNMT3A v009 Missense
S02 DNMT3A v010 Missense
S04 DNMT3A v011 Truncating
S08 DNMT3A v012 Missense
S09 DNMT3A v013 Splice
S03 RUNX1 v014 Truncating
S06 RUNX1 v015 Missense
S07 RUNX1 v016 Splice
S10 RUNX1 v017 Truncating
S01 FLT3 v018 Inframe
S05 FLT3 v019 Missense
S09 FLT3 v020 Inframe
S12 FLT3 v021 Missense
S02 NPM1 v022 Inframe
S04 NPM1 v023 Inframe
S07 NPM1 v024 Inframe
S08 NPM1 v025 Inframe
S06 TET2 v026 Truncating
S10 TET2 v027 Missense
S11 TET2 v028 Truncating
S03 IDH1 v029 Missense
S09 IDH1 v030 Missense
S05 ASXL1 v031 Truncating
S12 ASXL1 v032 Splice
;
run;


%MGXOncoPrint(samples=work.absent,genes=work.genes,
 callability=work.callability,variants=work.variants,outdir=&BOOK_OUT,
 template_file="oncoprint_template.sas");
data work.mgx_negative_status;length test $40;test='Missing dataset';status=&MGX_ONCOPRINT_RC;output;run;
data work.bad_variants;set work.variants;if _n_=1 then mutation_class='Unsupported';run;
%MGXOncoPrint(samples=work.samples,genes=work.genes,
 callability=work.callability,variants=work.bad_variants,outdir=&BOOK_OUT,
 template_file="oncoprint_template.sas");
data work._mgx_status;length test $40;test='Invalid mutation class';status=&MGX_ONCOPRINT_RC;run;
proc append base=work.mgx_negative_status data=work._mgx_status;run;
data work.bad_calls;set work.callability;if sample='S01' and gene_index=8 then callable=0;run;
%MGXOncoPrint(samples=work.samples,genes=work.genes,
 callability=work.bad_calls,variants=work.variants,outdir=&BOOK_OUT,
 template_file="oncoprint_template.sas");
data work._mgx_status;length test $40;test='Event in unassessable pair';status=&MGX_ONCOPRINT_RC;run;
proc append base=work.mgx_negative_status data=work._mgx_status;run;
