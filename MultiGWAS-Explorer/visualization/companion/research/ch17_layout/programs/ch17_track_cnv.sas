/* Real GRCh37 annotation; three explicitly selected CYP26B1 transcripts. */
data work.track_gtf;
  length type $10 genesymbol $12 ensembl_transcript $35;
  infile "track_annotation.tsv" dsd dlm='09'x firstobs=2 truncover;
  input chr st end type :$10. genesymbol :$12. ensembl_transcript :$35. protein_coding;
  span=.;
run;
%let local_gtf_ods_html5=1;
ods listing close;
ods html5 path="&BOOK_OUT" gpath="&BOOK_OUT" file='tracks_native.html'
 options(bitmap_mode='LINK') image_dpi=300 style=BookPublication;
/* Invented association examples, not observed CNV calls or PGC estimates. */
data work.cnv_input;
  length snp $12 variant_kind $3;
  input snp :$12. variant_kind :$3. left right z;
  chr=2;bp=mean(left,right);
  if variant_kind='CNV' then span=right-left;else span=.;
  p=2*sdf('normal',abs(z));
  datalines;
SNV_A SNV 72351200 72351200 -0.8
SNV_B SNV 72354500 72354500 1.6
CNV_A CNV 72357000 72361500 -4.4
SNV_C SNV 72359800 72359800 -3.1
SNV_D SNV 72363000 72363000 2.4
CNV_B CNV 72363500 72370500 4.9
SNV_E SNV 72367000 72367000 2.9
SNV_F SNV 72371200 72371200 -1.5
CNV_C CNV 72372500 72377000 -3.6
SNV_G SNV 72374500 72374500 -2.6
SNV_H SNV 72378500 72378500 0.7
;
run;
%SNP_Local_Manhattan_With_GTF(gwas_dsd=work.cnv_input,chr_var=chr,
 AssocPVars=p,ZscoreVars=z,SNP_IDs=2:72350000:72380000,
 SNP_Var=snp,Pos_Var=bp,Variant_Length_Var=span,
 dist2snp=0,gtf_dsd=work.track_gtf,focus_on_transcript=0,
 gwas_labels_in_order=Invented_association,
 design_width=760,design_height=400,image_dpi=300,
 dotsize=5,highlow_line_cmd=%str(thickness=4 color=darkorange pattern=solid),
 dist2sep_genes=50000,pct4neg_y=0.5,grp_font_size=9,
 fig_fmt=png,makedotheatmap=1,makeheatmapdotintooneline=0,
 heatmap_min_neg_val=-5,heatmap_max_pos_val=5,
 heatmap_Neg_rangealtcolormodel=CX2166AC CX92C5DE CXF7F7F7,
 heatmap_Pos_rangealtcolormodel=CXF7F7F7 CXF4A582 CXB2182B,
 heatmap_legend_title=%str(Invented effect-direction statistic),verbose=1);
data work.cnv_plot;set work.final;run;
options notes;
ods html5 close;
proc template;source BedGraph / file="&BOOK_OUT/cnv_intervals_native_template.sas";run;
filename tsrc "&BOOK_OUT/cnv_intervals_native_template.sas";
filename tout "&BOOK_OUT/cnv_intervals_template.sas";
data _null_;
 infile tsrc lrecl=32767;file tout lrecl=32767;length line $32767;input;line=_infile_;
 if index(lowcase(line),'entrytitle') then line='entrytitle "Invented SNV and CNV association examples; real CYP26B1 annotation";';
 put line;
run;
proc template;
%include tout;
run;
title;footnote;
ods listing gpath="&BOOK_OUT" image_dpi=300 style=BookPublication;
ods graphics / reset=index imagename='cnv_intervals' imagefmt=png width=7.6in height=4in;
proc sgrender data=work.final template=BedGraph;dynamic _chr='CHR';run;
ods graphics / reset=index imagename='cnv_intervals' imagefmt=svg width=7.6in height=4in;
proc sgrender data=work.final template=BedGraph;dynamic _chr='CHR';run;
ods listing close;

/* Presentation-only pass over the original macro-generated GTL. */
%book_pub_gtl(in=&BOOK_OUT/cnv_intervals_template.sas,out=&BOOK_OUT/publication_template.sas,chr=2,start=72350000,stop=72380000,step=5000,decimals=3);
title;footnote;
ods listing gpath="&BOOK_OUT" image_dpi=450 style=BookPublication;
ods graphics / reset=index imagename='cnv_intervals' imagefmt=png width=7.6in height=4.5in;
proc sgrender data=work.final template=BedGraph;dynamic _chr='CHR';run;
ods graphics / reset=index imagename='cnv_intervals' imagefmt=svg width=7.6in height=4.5in;
proc sgrender data=work.final template=BedGraph;dynamic _chr='CHR';run;
ods listing close;

/* Book-owned proportions: raw macro tables remain unchanged. */
%book_track_display(data=work.final,out=work.track_display,in=&BOOK_OUT/publication_template.sas,template=&BOOK_OUT/track_template.sas,factor=1.4088,body=14,single_marker=0,label_center_pct=,label_gap=0.293489957);
title;footnote;
ods listing gpath="&BOOK_OUT" image_dpi=450 style=BookPublication;
ods graphics / reset=index imagename='cnv_intervals' imagefmt=png width=7.6in height=4.5in;
proc sgrender data=work.track_display template=BedGraph;dynamic _chr='CHR';run;
/* Vector coordinates use native design resolution; avoid DPI-dependent SVG strokes. */
ods listing gpath="&BOOK_OUT" image_dpi=96 style=BookPublication;
ods graphics / reset=index imagename='cnv_intervals' imagefmt=svg width=7.6in height=4.5in;
proc sgrender data=work.track_display template=BedGraph;dynamic _chr='CHR';run;
ods listing close;

/* CNV display option 1: retain ST/END/POS; mask only the point layer. */
%BookCNVDisplay(data=work.track_display,out=work.cnv_midpoint_control,length_var=span,show_midpoints=1);
ods listing gpath="&BOOK_OUT" image_dpi=450 style=BookPublication;
ods graphics / reset=index imagename='cnv_midpoints_control' imagefmt=png width=7.6in height=4.5in;
proc sgrender data=work.cnv_midpoint_control template=BedGraph;dynamic _chr='CHR';run;
ods listing gpath="&BOOK_OUT" image_dpi=96 style=BookPublication;
ods graphics / reset=index imagename='cnv_midpoints_control' imagefmt=svg width=7.6in height=4.5in;
proc sgrender data=work.cnv_midpoint_control template=BedGraph;dynamic _chr='CHR';run;
ods listing close;

/* CNV display option 0: retain ST/END/POS; mask only the point layer. */
%BookCNVDisplay(data=work.track_display,out=work.cnv_display,length_var=span,show_midpoints=0);
ods listing gpath="&BOOK_OUT" image_dpi=450 style=BookPublication;
ods graphics / reset=index imagename='cnv_intervals' imagefmt=png width=7.6in height=4.5in;
proc sgrender data=work.cnv_display template=BedGraph;dynamic _chr='CHR';run;
ods listing gpath="&BOOK_OUT" image_dpi=96 style=BookPublication;
ods graphics / reset=index imagename='cnv_intervals' imagefmt=svg width=7.6in height=4.5in;
proc sgrender data=work.cnv_display template=BedGraph;dynamic _chr='CHR';run;
ods listing close;

/* Author-requested larger markers and tighter lower clearance. */
%BookChapter17Layout(data=work.cnv_display,out=work.chapter17_display,height=450,dots=7,bottom_pad=12,body=14);
/* Thicken only CNV intervals; preserve their ST/END and the masked point layer. */
data _null_;
 infile "&BOOK_OUT/ch17_template.sas" lrecl=32767 truncover;
 file "&BOOK_OUT/ch17_cnv_template.sas" lrecl=32767;
 length line $32767;input;line=_infile_;
 if prxmatch('/^highlowplot y=VAR4LOG10P high=END low=ST /i',line)
  then line=prxchange('s/thickness=[0-9.]+/thickness=8px/i',1,line);
 put line;
run;
proc template;%include "&BOOK_OUT/ch17_cnv_template.sas";run;
title;footnote;
ods listing gpath="&BOOK_OUT" image_dpi=450 style=BookPublication;
ods graphics / reset=index imagename='cnv_intervals' imagefmt=png width=7.6in height=4.5in;
proc sgrender data=work.chapter17_display template=BedGraph;dynamic _chr='CHR';run;
ods listing gpath="&BOOK_OUT" image_dpi=96 style=BookPublication;
ods graphics / reset=index imagename='cnv_intervals' imagefmt=svg width=7.6in height=4.5in;
proc sgrender data=work.chapter17_display template=BedGraph;dynamic _chr='CHR';run;
ods listing close;
