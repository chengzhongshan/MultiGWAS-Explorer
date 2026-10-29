/* Book-owned presentation controls. Original research macros are unchanged. */
proc template;
 define style styles.BookPublication;
 parent=styles.journal;
 style GraphFonts from GraphFonts /
  'GraphTitleFont'=("Arial",9pt)
  'GraphLabelFont'=("Arial",9pt)
  'GraphValueFont'=("Arial",8.5pt)
  'GraphDataFont'=("Arial",8.5pt)
  'GraphUnicodeFont'=("Arial",8.5pt);
 end;
run;

%macro book_pub_gtl(in=,out=,chr=,start=,stop=,step=,decimals=1);
%local tick ticks;
%let ticks=;
%if %length(&chr) %then %do;
 proc format;
 value bookmb
 %do tick=&start %to &stop %by &step;
  &tick="%sysfunc(strip(%sysfunc(putn(%sysevalf(&tick/1000000),12.&decimals))))"
  %let ticks=&ticks &tick;
 %end;
 ;
 run;
%end;
data _null_;
 infile "&in" lrecl=32767 truncover;
 file "&out" lrecl=32767;
 length line $32767;
 input;line=_infile_;
 /* Captions identify the panel; avoid repeating a long headline in the image. */
 if index(lowcase(line),'entrytitle') then delete;
 line=tranwrd(line,'color=lightgray pattern=1 thickness=5',
                   'color=CXCCCCCC pattern=1 thickness=0.5');
 %if %length(&chr) %then %do;
 if index(lowcase(line),'xaxisopts=(') then do;
  line=prxchange('s/xaxisopts=\(label="[^"]*"/xaxisopts=(/',1,line);
  line=tranwrd(line,'xaxisopts=(',
    'xaxisopts=(label="Chromosome &chr position (Mb)" ');
 end;
 line=tranwrd(line,'tickvalueformat=best32.',
   'tickvalueformat=bookmb. tickvaluelist=(&ticks)');
 %end;
 put line;
run;
proc template;%include "&out";run;
%mend;
