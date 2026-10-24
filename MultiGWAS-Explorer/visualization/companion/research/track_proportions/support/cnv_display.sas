/* Separate point coordinates from interval coordinates in original GTL. */
%macro BookCNVDisplay(data=work.track_display,out=work.cnv_display,
 length_var=span,show_midpoints=0);
%local book_original_columns book_cnv_compare;
proc contents data=&data out=work.book_cnv_columns noprint;run;
proc sql noprint;
 select name into :book_original_columns separated by ' '
 from work.book_cnv_columns order by varnum;
quit;
data &out;
 set &data;
 book_scatter_pos=pos;
 /* >1 matches the original macro's criterion for drawing HIGHLOW intervals. */
 if &show_midpoints=0 and grp='GWAS_Assoc_Signal' and not missing(old_y)
    and not missing(&length_var) and &length_var>1 then book_scatter_pos=.;
run;
proc compare base=&data compare=&out(keep=&book_original_columns) criterion=0 noprint;run;
%let book_cnv_compare=&sysinfo;
data work.cnv_display_checks;
 compare_sysinfo=&book_cnv_compare;show_midpoints=&show_midpoints;
 if compare_sysinfo ne 0 then put 'ERROR: CNV display changed original columns.';
run;
proc template;source BedGraph / file="&BOOK_OUT/cnv_resolved_template.sas";run;
data _null_;
 infile "&BOOK_OUT/cnv_resolved_template.sas" lrecl=32767 truncover;
 file "&BOOK_OUT/cnv_display_template.sas" lrecl=32767;
 length line $32767;input;line=_infile_;
 line=tranwrd(line,'scatterplot x=POS y=VAR4LOG10P /',
                   'scatterplot x=BOOK_SCATTER_POS y=VAR4LOG10P /');
 put line;
run;
proc template;%include "&BOOK_OUT/cnv_display_template.sas";run;
%mend;
