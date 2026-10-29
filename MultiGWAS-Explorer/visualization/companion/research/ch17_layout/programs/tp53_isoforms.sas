/* Real UniProt protein isoforms, aligned to P04637-1 reference coordinates.
   Canonical GDC sites are annotation guides, not isoform-specific tumor calls. */
data work.isoform_features;
 length type $10 genesymbol $20;
 infile "tp53_isoform_features.tsv" dsd dlm='09'x firstobs=2 truncover;
 input chr st end type :$10. genesymbol :$20. protein_coding;
run;
data work.isoform_metadata;
 length isoform $20 isoform_name $30 seqhash $64;
 infile "tp53_isoforms.tsv" dsd dlm='09'x firstobs=2 truncover;
 input isoform :$20. isoform_name :$30. protein_length canonical_start canonical_end
       alternative_start alternative_end seqhash :$64.;
 isoform_order=_n_;book_alt_y=-isoform_order;
 shared_end=canonical_end;
 if not missing(alternative_start) then shared_end=alternative_start-1;
 book_label_x=canonical_start;book_label_y=-isoform_order+0.28;
 length book_iso_label $90;
 book_iso_label=catx(' ',isoform_name,cats('(',isoform,','),cats(protein_length),'aa)');
run;
data work.mutation_long;
 length marker $12 panel $8 mutation_type $16;
 infile "protein_occurrences.tsv" dsd dlm='09'x firstobs=2 truncover;
 input marker :$12. chr pos panel :$8. mutation_type :$16.;
 if pos in (175,220,248,273);present=1;se=.;p=.5;
run;
%let local_gtf_ods_html5=1;
ods listing close;
ods html5 path="&BOOK_OUT" gpath="&BOOK_OUT" file='isoform_macro.html'
 options(bitmap_mode='LINK') image_dpi=96 style=BookPublication;
%BookGDCMutationProteinTrack(longformat_muts=work.mutation_long,gwas_dsd_var=panel,
 query_rsid_or_genesymbol=1:1:393,dist2snp_or_gene=0,
 subset_gwas_dsd=work.mutation_wide,gwas_chr_var=chr,
 gwas_p_var=p,gwas_beta_var=present,gwas_se_var=se,
 gwas_snp_var=marker,gwas_pos_var=pos,GenerateGeneManhattanPlot=1,
 SNPs2label_scatterplot_dots=R175H Y220C R248Q R273H,
 gtf_dsd=work.isoform_features,design_width=760,design_height=620,
 barthickness=14,dotsize=7,dist2sep_genes=10,
 yoffset4max_drawmarkersontop=0.15,Yoffset4textlabels=1.1,
 text_rotate_angle=90,auto_rotate2zero=0,adj_spaces_among_top_snps=0,
 gwas_labels_in_order=,outfigfmt=png);
options notes;ods html5 close;
data work.isoform_native_plot;set work.final;run;
/* Retain the macro's annotation x coordinates; place each isoform on its own row. */
proc contents data=work.isoform_native_plot out=work.iso_columns noprint;run;
proc sql noprint;
 select name into :iso_yvars separated by ' ' from work.iso_columns
 where prxmatch('/^VAR4LOG10P[0-9]+$/i',strip(name)) order by varnum;
 select name into :iso_groups separated by ' ' from work.iso_columns
 where prxmatch('/^GRP[0-9]+$/i',strip(name)) order by varnum;
 select name into :iso_allvars separated by ' ' from work.iso_columns order by varnum;
quit;
data work.isoform_display;
 set work.isoform_native_plot;
 array gy {*} &iso_yvars;array gg {*} &iso_groups;
 array book_iso_y {*} book_iso_y1-book_iso_y%sysfunc(countw(&iso_yvars));
 do i=1 to dim(gy);
  if not missing(gy[i]) and gy[i]<0 and not missing(gg[i])
   then book_iso_y[i]=-input(scan(gg[i],2,'-'),best32.);
 end;
 drop i;
run;
proc compare base=work.isoform_native_plot
 compare=work.isoform_display(keep=&iso_allvars) criterion=0 noprint;run;
data work.isoform_compare;compare_sysinfo=&sysinfo;run;
data work.isoform_canvas;set work.isoform_display work.isoform_metadata;run;
proc template;source BedGraph / file="&BOOK_OUT/isoform_original_template.sas";run;
data _null_;
 infile "&BOOK_OUT/isoform_original_template.sas" lrecl=32767 truncover;
 file "&BOOK_OUT/isoform_final_template.sas" lrecl=32767;
 length line $32767;retain skip_axis 0;input;line=_infile_;
 if index(lowcase(line),'layout lattice /') then do;
  put 'legenditem type=line name="shared" / label="Sequence identical to p53alpha" lineattrs=(color=CX0072B2 thickness=8);';
  put 'legenditem type=line name="alternate" / label="Alternative C-terminal sequence" lineattrs=(color=CXD55E00 thickness=8);';
  put 'legenditem type=line name="site" / label="Projected canonical GDC site" lineattrs=(color=CX888888 pattern=shortdash);';
 end;
 if index(lowcase(line),'layout overlay / yaxisopts=') then do;
  put 'layout overlay / yaxisopts=(display=none offsetmin=0.01 offsetmax=0.01 linearopts=(viewmin=-9.3 viewmax=0.3))';
  put 'xaxisopts=(label="Aligned position in canonical p53alpha (P04637-1)" offsetmin=0.015 offsetmax=0.015 linearopts=(viewmin=1 viewmax=393 tickvaluelist=(1 50 100 150 200 250 300 350 393)));';
  put 'referenceline x=175 / lineattrs=(color=CX888888 pattern=shortdash);';
  put 'referenceline x=220 / lineattrs=(color=CX888888 pattern=shortdash);';
  put 'referenceline x=248 / lineattrs=(color=CX888888 pattern=shortdash);';
  put 'referenceline x=273 / lineattrs=(color=CX888888 pattern=shortdash);';
  put 'drawtext textattrs=(size=9) "R175H" / x=175 y=94 xspace=datavalue yspace=wallpercent anchor=center rotate=90;';
  put 'drawtext textattrs=(size=9) "Y220C" / x=220 y=94 xspace=datavalue yspace=wallpercent anchor=center rotate=90;';
  put 'drawtext textattrs=(size=9) "R248Q" / x=248 y=94 xspace=datavalue yspace=wallpercent anchor=center rotate=90;';
  put 'drawtext textattrs=(size=9) "R273H" / x=273 y=94 xspace=datavalue yspace=wallpercent anchor=center rotate=90;';
  /* A full-length exon-like interval can be collapsed to fake coordinates by
     the inherited annotation transformer. Draw both sequence layers from the
     independently checked UniProt interval table; retain native rows for audit. */
  put 'highlowplot y=BOOK_ALT_Y low=CANONICAL_START high=SHARED_END / lineattrs=(color=CX0072B2 thickness=14px);';
  put 'highlowplot y=BOOK_ALT_Y low=ALTERNATIVE_START high=ALTERNATIVE_END / lineattrs=(color=CXD55E00 thickness=14px);';
  put 'scatterplot x=BOOK_LABEL_X y=BOOK_LABEL_Y / markercharacter=BOOK_ISO_LABEL markercharacterposition=topright markercharacterattrs=(size=9 color=black);';
  skip_axis=1;delete;
 end;
 if skip_axis then do;if index(lowcase(line),'offsetmin=') then skip_axis=0;delete;end;
 if index(lowcase(line),'seriesplot x=') then delete;
 if index(lowcase(line),'scatterplot x=') then delete;
 if prxmatch('/(continuouslegend|textplot|referenceline y=|entrytitle|drawline)/i',line) then delete;
 if index(lowcase(line),'sidebar / align=bottom;') then do;
  put line;put 'discretelegend "shared" "alternate" "site" / border=false valueattrs=(size=8);';delete;
 end;
 put line;
run;
proc template;%include "&BOOK_OUT/isoform_final_template.sas";run;
title;footnote;
ods listing gpath="&BOOK_OUT" image_dpi=450 style=BookPublication;
ods graphics / reset=index imagename='tp53_isoforms' imagefmt=png width=7.6in height=6.2in;
proc sgrender data=work.isoform_canvas template=BedGraph;run;
ods listing gpath="&BOOK_OUT" image_dpi=96 style=BookPublication;
ods graphics / reset=index imagename='tp53_isoforms' imagefmt=svg width=7.6in height=6.2in;
proc sgrender data=work.isoform_canvas template=BedGraph;run;
ods listing close;
