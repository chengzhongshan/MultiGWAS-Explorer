%macro rank4grps(grps,dsdout);
data &dsdout;
length grps $300.;
%do i=1 %to %sysfunc(countw(%superq(grps),%str( )));
  %let gval=%qscan(%superq(grps),&i,%str( ));
  grps="&gval";num_grps=&i;
  output;
%end;
run;
proc sort data=&dsdout;by grps;
data &dsdout;
set &dsdout;
char_ord=_n_;
run;
proc sort data=&dsdout;by num_grps;
run;

%mend;
