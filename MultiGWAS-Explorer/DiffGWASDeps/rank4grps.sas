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

/*Demo:

%rank4grps(
grps=rs8116534 rs472481 rs555336963 rs148143613 rs2924725 rs5927942,
dsdout=z
);
proc print;run;
*Note: the new variable char_ord is a index for extracting these grps in the alphabet sorted array of these grps;
*This means when the grps is sorted with the default mode by alphabet order, the char_ord is an index to get back the original order of grps;
*based on the sorted grps in alphabet order;

*/
