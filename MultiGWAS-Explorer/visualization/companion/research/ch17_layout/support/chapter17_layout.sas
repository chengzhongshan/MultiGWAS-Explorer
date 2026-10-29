/* Scoped display pass. Original macros, data coordinates and palettes remain intact. */
%macro BookChapter17Layout(data=,out=work.chapter17_display,height=400,dots=6,
 bottom_pad=12,body=14);
%local annvars allvars lower_seen compare_status;
proc template;source BedGraph / file="&BOOK_OUT/ch17_before_template.sas";run;
proc contents data=&data out=work.ch17_columns noprint;run;
proc sql noprint;
 select name into :annvars separated by ' ' from work.ch17_columns
 where prxmatch('/^BOOK_ANN[0-9]+$/i',strip(name)) order by varnum;
 select name into :allvars separated by ' ' from work.ch17_columns order by varnum;
quit;
data &out;set &data;run;
proc compare base=&data compare=&out(keep=&allvars) criterion=0 noprint;run;
%let compare_status=&sysinfo;
data work.chapter17_compare;compare_sysinfo=&compare_status;run;
%if &compare_status ne 0 %then %do;
 %put ERROR: Chapter 17 display copy changed original values.;
 %return;
%end;
/* Find the lowest annotation center, not a generic negative-axis fraction. */
data _null_;
 set &data end=last;
 array a {*} &annvars;
 retain lowest 0;
 do i=1 to dim(a);if not missing(a[i]) then lowest=min(lowest,a[i]);end;
 if last then call symputx('lower_seen',lowest,'L');
run;
/* Preserve the positive scale. Reserve a small clearance below thick body bars. */
data work.chapter17_layout_parameters;
 infile "&BOOK_OUT/ch17_before_template.sas" lrecl=32767 truncover end=last;
 length line $32767;retain upper offset_top;
 input;line=lowcase(_infile_);
 rx=prxparse('/viewmax=([-0-9.]+) viewmin=(-[0-9.]+)/');
 if prxmatch(rx,line) then do;
  call prxposn(rx,1,s,n);upper=input(substr(line,s,n),best32.);
 end;
 rx2=prxparse('/offsetmin=[0-9.]+ offsetmax=([0-9.]+)/');
 if prxmatch(rx2,line) then do;
  call prxposn(rx2,1,s,n);offset_top=input(substr(line,s,n),best32.);
 end;
 if last then do;
  lowest=&lower_seen;height=&height;marker_size=&dots;bottom_padding=&bottom_pad;
  usable=(height-70)*(1-offset_top-0.01);
  pad=max(bottom_padding,&body/2+4);
  lower=lowest-pad*(upper-lowest)/(usable-pad);
  if nmiss(upper,offset_top,lower) or usable<=pad then do;
   put 'ERROR: Unsupported chapter 17 GTL axis layout.';stop;
  end;
  call symputx('ch17_lower',lower,'L');output;
 end;
 keep upper lower lowest height marker_size bottom_padding;
run;
data _null_;
 infile "&BOOK_OUT/ch17_before_template.sas" lrecl=32767 truncover;
 file "&BOOK_OUT/ch17_template.sas" lrecl=32767;
 length line $32767;input;line=_infile_;
 line=prxchange("s/designheight=[0-9.]+/designheight=&height/i",1,line);
 if index(lowcase(line),'yaxisopts=') then do;
  line=prxchange("s/viewmin=-[0-9.]+/viewmin=&ch17_lower/i",1,line);
  line=prxchange('s/offsetmin=[0-9.]+/offsetmin=0.01/i',1,line);
 end;
 /* Only the actual association/occurrence point layer gets larger markers. */
 if prxmatch('/scatterplot x=(POS|BOOK_SCATTER_POS|BOOK_MARKER_POS) y=VAR4LOG10P /i',line)
  then line=prxchange("s/size=[0-9.]+(?:px|pt)?/size=&dots.px/i",1,line);
 put line;
run;
proc template;%include "&BOOK_OUT/ch17_template.sas";run;
%mend;
