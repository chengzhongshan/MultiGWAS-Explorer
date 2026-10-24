/* Book-owned display pass: preserve every original macro-generated column. */
%macro book_track_display(data=,out=work.track_display,in=,template=,
 factor=0.5,body=16,single_marker=0,label_center_pct=,label_gap=0);
%local annvars allvars nann i compare_status;
/* SOURCE materializes the compiled template, including resolved local macro
   variables from book_pub_gtl (chromosome and tick list). */
proc template;source BedGraph / file="&BOOK_OUT/track_resolved_template.sas";run;
%let in=&BOOK_OUT/track_resolved_template.sas;
proc contents data=&data out=work.book_columns noprint;run;
proc sql noprint;
 select name into :annvars separated by ' ' from work.book_columns
 where prxmatch('/^VAR4LOG10P\d+$/i',strip(name)) order by varnum;
 select name into :allvars separated by ' ' from work.book_columns order by varnum;
quit;
%let nann=%sysfunc(countw(&annvars));
%if &single_marker %then %do;
 proc sql;
  create table work.book_group_order as select name from dictionary.columns
   where libname='WORK' and memname='MUTATION_WIDE' and upcase(substr(name,1,2))='P_'
   order by varnum;
 quit;
 data work.book_group_order;set work.book_group_order;pi=_n_;run;
 proc sql;
  create table work.book_occurrence_key as select distinct a.pos,b.pi
   from work.mutation_long a inner join work.book_group_order b
   on lowcase(a.panel)=lowcase(substr(b.name,3));
 quit;
%end;
data &out;
 %if &single_marker %then %do;
 if _n_=1 then do;
  declare hash occurrence(dataset:'work.book_occurrence_key');
  occurrence.defineKey('pi','pos');occurrence.defineDone();
 end;
 %end;
 set &data;
 array original_ann {*} &annvars;
 array display_ann {*} book_ann1-book_ann&nann;
 do book_i=1 to dim(original_ann);
  display_ann[book_i]=original_ann[book_i];
  if not missing(original_ann[book_i]) and original_ann[book_i]<0 then
   display_ann[book_i]=&factor*original_ann[book_i];
 end;
 book_gene_y=_y_;
 if not missing(_y_) and _y_<0 then book_gene_y=&factor*_y_;
 if not missing(book_gene_y) and book_gene_y<0 and grp_label ne ''
    and not missing(display_ann[1]) and display_ann[1]<0 then
  book_gene_y=max(book_gene_y,display_ann[1]+&label_gap);
 %if &single_marker %then %do;
 book_marker_pos=.;
 if grp='GWAS_Assoc_Signal' and not missing(old_y) and occurrence.check()=0 then book_marker_pos=pos;
 %end;
 drop book_i;
run;
/* Exact check of all original columns; extra display columns are excluded. */
proc compare base=&data compare=&out(keep=&allvars) criterion=0 noprint;run;
%let compare_status=&sysinfo;
data work.track_checks;
 compare_sysinfo=&compare_status;annotation_scale=&factor;body_thickness=&body;
 single_marker=&single_marker;
 if compare_sysinfo ne 0 then put 'ERROR: Original plotting columns changed.';
run;
%if %length(&label_center_pct) %then %do;
 proc sort data=&data(keep=_target_snp_ newpos where=(_target_snp_ ne ''))
  out=work.book_labels nodupkey;by _target_snp_;run;
%end;
data _null_;
 infile "&in" lrecl=32767 truncover;file "&template" lrecl=32767;
 length line token $32767;retain rx;
 if _n_=1 then rx=prxparse('/viewmin=(-[0-9.]+)/i');
 input;line=_infile_;
 /* Compress the negative annotation coordinates only, on the same x axis. */
 if prxmatch(rx,line) then do;
  call prxposn(rx,1,book_start,book_length);
  token=substr(line,book_start,book_length);
  line=tranwrd(line,cats('viewmin=',token),cats('viewmin=',put(input(token,best32.)*&factor,best32.)));
 end;
 if index(lowcase(line),'seriesplot x=') then do;
  %do i=1 %to &nann;
  line=tranwrd(line,"y=%upcase(%scan(&annvars,&i)) /","y=BOOK_ANN&i /");
  %end;
  line=prxchange("s/thickness=[0-9.]+(?:px|pt|in|cm|mm)?/thickness=&body.px/",1,line);
 end;
 line=tranwrd(line,'y=_Y_ /','y=BOOK_GENE_Y /');
 %if &single_marker %then %do;
 line=tranwrd(line,'scatterplot x=POS y=VAR4LOG10P /','scatterplot x=BOOK_MARKER_POS y=VAR4LOG10P /');
 %end;
 %if %length(&label_center_pct) %then %do;
 if index(lowcase(line),'textplot x=newpos') then do;
  /* Center the complete rotated string inside the reserved wall band. */
  do book_label_i=1 to book_nlabels;
   set work.book_labels point=book_label_i nobs=book_nlabels;
   put 'drawtext textattrs=(size=8) "' _target_snp_ +(-1)
    '" / x=' newpos best32. " y=&label_center_pct xspace=datavalue yspace=wallpercent anchor=center rotate=90;";
  end;
  delete;
 end;
 %end;
 put line;
run;
proc template;%include "&template";run;
%mend;
