/* Presence-only records. P=.5 is an internal drawing placeholder; no test. */
%let local_gtf_ods_html5=1;
ods listing close;
ods html5 path="&BOOK_OUT" gpath="&BOOK_OUT" file='mutation_native.html'
 options(bitmap_mode='LINK') image_dpi=300 style=BookPublication;
/* UniProt P04637: functional regions use the wrapper's exon drawing type. */
data work.protein_features;
 length genesymbol $20 type $10;
 infile "tp53_protein_features.tsv" dsd dlm='09'x firstobs=2 truncover;
 input chr st end type :$10. genesymbol :$20. protein_coding;
run;

/* Source: frozen open GDC somatic variants. P=.5 is a drawing placeholder. */
data work.mutation_long;
 length marker $12 panel $8 mutation_type $16;
 infile "protein_occurrences.tsv" dsd dlm='09'x firstobs=2 truncover;
 input marker :$12. chr pos panel :$8. mutation_type :$16.;
 present=1;se=.;p=.5;
run;
%BookGDCMutationProteinTrack(longformat_muts=work.mutation_long,gwas_dsd_var=panel,
 query_rsid_or_genesymbol=1:1:393,dist2snp_or_gene=0,
 subset_gwas_dsd=work.mutation_wide,gwas_chr_var=chr,
 gwas_p_var=p,gwas_beta_var=present,gwas_se_var=se,
 gwas_snp_var=marker,gwas_pos_var=pos,GenerateGeneManhattanPlot=1,
 SNPs2label_scatterplot_dots=W53Ter R175H Y220C R248Q R273H R342Ter,
 gtf_dsd=work.protein_features,design_width=760,design_height=420,
 barthickness=10,dotsize=7,dist2sep_genes=10,
 yoffset4max_drawmarkersontop=0.30,Yoffset4textlabels=1.1,
 text_rotate_angle=90,auto_rotate2zero=0,adj_spaces_among_top_snps=0,
 gwas_labels_in_order=,outfigfmt=png);
options notes;
ods html5 close;
/* Add categorical colors to the original drawing records without changing sites. */
proc sort data=work.mutation_long(keep=pos mutation_type)
          out=work.mutation_type_key nodupkey;
  by pos;
run;
data work.mutation_native_plot;
  set work.final;
  if grp='GWAS_Assoc_Signal' and not missing(old_y) then AssocGrp=old_y;
run;
proc sql;
  create table work.mutation_plot as
  select p.*, k.mutation_type
  from work.mutation_native_plot as p
  left join work.mutation_type_key as k
    on p.grp='GWAS_Assoc_Signal' and not missing(p.old_y)
       and p.pos between k.pos and k.pos+1;
quit;

data work.mutation_plot;set work.mutation_plot;label pos='Amino-acid position in TP53 (P04637)' newpos='Amino-acid position in TP53 (P04637)';run;
proc template;source BedGraph / file="&BOOK_OUT/gdc_tp53_protein_native_template.sas";run;
filename tsrc "&BOOK_OUT/gdc_tp53_protein_native_template.sas";
filename tout "&BOOK_OUT/gdc_tp53_protein_template.sas";
data _null_;infile tsrc lrecl=32767;file tout lrecl=32767;length line $32767;input;line=_infile_;
 if index(lowcase(line),'entrytitle') then line='entrytitle "GDC TP53 somatic variants: TCGA-BRCA and TCGA-LUAD";';
 if index(lowcase(line),'continuouslegend "sc"') then delete;
 line=tranwrd(line,'label="-log10(P)"','label="Project presence"');
 if index(lowcase(line),'discretelegend "sc"') then delete;
 line=tranwrd(line,'display=(tickvalues label)','display=(label)');
 if index(lowcase(line),'layout lattice /') then do;
   put 'discreteattrmap name="MutationTypes" / ignorecase=true;';
   put 'value "Missense" / markerattrs=(color=CX0072B2 symbol=circlefilled);';
   put 'value "Nonsense" / markerattrs=(color=CXD55E00 symbol=circlefilled);';
   put 'value "Frameshift" / markerattrs=(color=CXCC79A7 symbol=circlefilled);';
   put 'value "Splice site" / markerattrs=(color=CX009E73 symbol=circlefilled);';
   put 'enddiscreteattrmap;';
   put 'discreteattrvar attrvar=MutationTypeColor var=mutation_type attrmap="MutationTypes";';
 end;
 line=tranwrd(line,'markercolorgradient=heatmap_var_attrvar','group=MutationTypeColor');
 line=tranwrd(line,'group=AssocGrp name="sc"','group=MutationTypeColor name="sc"');
 if index(lowcase(line),'scatterplot x=pos y=var4log10p') then
   line=tranwrd(line,' / ',' / includemissinggroup=false ');
 if index(lowcase(line),'endsidebar;') then
   put 'discretelegend "sc" / title="GDC consequence" across=2 border=false;';
 line=tranwrd(line,'xaxisopts=(', 'xaxisopts=(label="Amino-acid position in TP53 (P04637)" ');
 if index(lowcase(line),'layout overlay / yaxisopts=') then
   line=prxchange('s/offsetmax=[0-9.]+/offsetmax=0.30/',1,line);
 if index(lowcase(line),'textplot x=newpos') then
   line=prxchange('s/POSITIONOFFSETY=[0-9.]+/POSITIONOFFSETY=1.1/',1,line);
 put line;run;
proc template;%include tout;run;
title;footnote;
ods listing gpath="&BOOK_OUT" image_dpi=300 style=BookPublication;
ods graphics / reset=index imagename='gdc_tp53_protein' imagefmt=png width=7.6in height=4.2in;
proc sgrender data=work.mutation_plot template=BedGraph;dynamic _chr='CHR';run;
ods graphics / reset=index imagename='gdc_tp53_protein' imagefmt=svg width=7.6in height=4.2in;
proc sgrender data=work.mutation_plot template=BedGraph;dynamic _chr='CHR';run;
ods listing close;

/* Presentation-only pass over the original macro-generated GTL. */
%book_pub_gtl(in=&BOOK_OUT/gdc_tp53_protein_template.sas,out=&BOOK_OUT/publication_template.sas,chr=,start=0,stop=0,step=1,decimals=0);
title;footnote;
ods listing gpath="&BOOK_OUT" image_dpi=450 style=BookPublication;
ods graphics / reset=index imagename='gdc_tp53_protein' imagefmt=png width=7.6in height=4.2in;
proc sgrender data=work.mutation_plot template=BedGraph;dynamic _chr='CHR';run;
ods graphics / reset=index imagename='gdc_tp53_protein' imagefmt=svg width=7.6in height=4.2in;
proc sgrender data=work.mutation_plot template=BedGraph;dynamic _chr='CHR';run;
ods listing close;

/* Book-owned proportions: raw macro tables remain unchanged. */
%book_track_display(data=work.mutation_plot,out=work.track_display,in=&BOOK_OUT/publication_template.sas,template=&BOOK_OUT/track_template.sas,factor=0.45,body=16,single_marker=1,label_center_pct=84.9,label_gap=0.746876254);
title;footnote;
ods listing gpath="&BOOK_OUT" image_dpi=450 style=BookPublication;
ods graphics / reset=index imagename='gdc_tp53_protein' imagefmt=png width=7.6in height=4.2in;
proc sgrender data=work.track_display template=BedGraph;dynamic _chr='CHR';run;
/* Vector coordinates use native design resolution; avoid DPI-dependent SVG strokes. */
ods listing gpath="&BOOK_OUT" image_dpi=96 style=BookPublication;
ods graphics / reset=index imagename='gdc_tp53_protein' imagefmt=svg width=7.6in height=4.2in;
proc sgrender data=work.track_display template=BedGraph;dynamic _chr='CHR';run;
ods listing close;

/* Compact occurrence rows and larger markers; no new mutation records. */
%BookMutationTrackDisplay(data=work.track_display,out=work.mutation_compact_display,height=300,dots=11,label_lift=0);
ods listing gpath="&BOOK_OUT" image_dpi=450 style=BookPublication;
ods graphics / reset=index imagename='gdc_tp53_protein' imagefmt=png width=7.6in height=3in;
proc sgrender data=work.mutation_compact_display template=BedGraph;dynamic _chr='CHR';run;
ods listing gpath="&BOOK_OUT" image_dpi=96 style=BookPublication;
ods graphics / reset=index imagename='gdc_tp53_protein' imagefmt=svg width=7.6in height=3in;
proc sgrender data=work.mutation_compact_display template=BedGraph;dynamic _chr='CHR';run;
ods listing close;

/* Author-requested larger markers and tighter lower clearance. */
%BookChapter17Layout(data=work.mutation_compact_display,out=work.chapter17_display,height=300,dots=14,bottom_pad=20,body=16);
title;footnote;
ods listing gpath="&BOOK_OUT" image_dpi=450 style=BookPublication;
ods graphics / reset=index imagename='gdc_tp53_protein' imagefmt=png width=7.6in height=3in;
proc sgrender data=work.chapter17_display template=BedGraph;dynamic _chr='CHR';run;
ods listing gpath="&BOOK_OUT" image_dpi=96 style=BookPublication;
ods graphics / reset=index imagename='gdc_tp53_protein' imagefmt=svg width=7.6in height=3in;
proc sgrender data=work.chapter17_display template=BedGraph;dynamic _chr='CHR';run;
ods listing close;
