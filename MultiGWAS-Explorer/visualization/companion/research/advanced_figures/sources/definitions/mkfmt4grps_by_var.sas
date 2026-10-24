%macro mkfmt4grps_by_var(
grpdsd,
grp_var,
by_var,/*Empty value is possible for this var, which is handy 
when there would be duplciates of grp_var after sorting by 
grp_var and by_var for the by_var to make a new format*/
outfmt4numgrps,
outfmt4chargrps,
dsd4fmt=dsd4fmt, /*This dataset can be used to get fmt info*/
numstarter=1 /*This will restrict the starter as 1 or other supplied number;
such as 0, which would be helpful when it is necessary to use different range
for making the format*/
);
%local rnd;

proc sort data=&grpdsd(keep=&grp_var &by_var) nodupkeys out=grpdsd;
*Make sure to sort it first by &by_var, and then by &grp_var;
*This was wrong, as some grps will be excluded by nodupkeys sorting with &by_var first;
by &grp_var &by_var;
run;

proc sql noprint;
select max(ndup_grps) into: maxdups
from 
(select count(*) as ndup_grps
from grpdsd
group by &grp_var)
;

%if &maxdups >1 %then %do;
		   %put Please evaluate the data set grpdsd for the combination between two vars, including &grp_var and &by_var;
		   %put It is necessary to have no duplicates in the first grp_var: &grp_var;
		   %put Currently there are duplicates for grp_var &grp_var after running proc sort nodupkeys by the two vars;
		   %put You may try to replace the two vars by using grp_var and grp_var, which means to use grp_var to replace the by_var;
		   %abort 255;
%end;

*The above dataset may still contains duplicates of &grp_var;
*Need to keep all unique &by_var and &grp_var;
/* proc sort data=grpdsd nodupkeys;by &grp_var; */
/* proc sort data=grpdsd;by &by_var &grp_var; */
/* run; */

*Need to sort grps by the by_var for making formats;
*Only when the by_var is not empty, run this codes;
%if %length(&by_var)>0 %then %do;
proc sort data=grpdsd;by &by_var;run;
%end;

/* Use a timestamp-derived suffix so this helper does not depend on
   a separately compiled random-number macro. */
%let rnd=%sysfunc(floor(%sysfunc(datetime())));
%let rnd=&rnd._;
data fmt4numgrps&rnd;
set grpdsd;
fmtname = "&outfmt4numgrps";
type = "I";
label = _N_ +  &numstarter - 1;
rename &grp_var = start;
run;
proc format cntlin=fmt4numgrps&rnd;
run;
quit;
data fmt4chargrps&rnd;
set grpdsd(drop=&by_var);
*Do not use &by_var here to asign group specific order for start;
*This will be a bug when there are duplicates in the &by_var;
/* by &by_var; */
start = _N_ + &numstarter - 1;
rename &grp_var=label;
fmtname = "&outfmt4chargrps";
type = "n";
run;
proc format cntlin=fmt4chargrps&rnd;
run;
%if %length(&dsd4fmt)>0 %then %do;
data &dsd4fmt;
set fmt4chargrps&rnd;
run;
%end;

%mend;
