%macro get_top_signal_within_dist(
dsdin=,
grp_var=,
signal_var=,
select_smallest_signal=1,
pos_var=,
pos_dist_thrshd=, /*Total exclusion span in bp; the macro keeps one lead signal within +/- 0.5*pos_dist_thrshd*/
dsdout=,
signal_thrshd=1 /*filter the input dsdin by association P, i.e, &signal_var <= &signal_thrshd*/
);

/*
Greedy in-group top-signal selection with a distance exclusion window.

Compared with the older self-join / SQL approach, this version is much faster
when many loci pass the threshold because it:

1. filters candidates first
2. sorts once by group and signal strength
3. keeps selected intervals for the current group only
4. greedily accepts the next best non-overlapping signal

The exclusion window matches the earlier behavior:
  dis_st = pos - pos_dist_thrshd * 0.5
  dis_end = pos + pos_dist_thrshd * 0.5

So a value such as 1e6 means "keep one lead signal per 1 Mb span"
rather than "use a 1 Mb half-window on each side".
*/

data _top_signal_candidates;
length Key $200.;
set &dsdin;
where not missing(&signal_var) and &signal_var > 0 and &signal_var <= &signal_thrshd;
dis_st=&pos_var-&pos_dist_thrshd*0.5;
dis_end=&pos_var+&pos_dist_thrshd*0.5;
Key=catx(':',&grp_var,&pos_var);
run;

%put Warning: duplicate records by Key and &signal_var will be kept for only one record.;
proc sort data=_top_signal_candidates dupout=_dup_records_ out=_top_signal_candidates nodupkeys;
%if &select_smallest_signal=1 %then %do;
  by &grp_var &pos_var &signal_var;
%end;
%else %do;
  by &grp_var &pos_var descending &signal_var;
%end;
run;
%put Warning: duplicate records by Key and &signal_var are saved into the data set _dup_records_.;

proc sort data=_top_signal_candidates;
%if &select_smallest_signal=1 %then %do;
  by &grp_var &signal_var &pos_var;
%end;
%else %do;
  by &grp_var descending &signal_var &pos_var;
%end;
run;

data &dsdout;
  array _sel_dis_st[50000] _temporary_;
  array _sel_dis_end[50000] _temporary_;
  retain _sel_n 0;
  length _keep 8 _i 8;

  set _top_signal_candidates;
  by &grp_var;

  if first.&grp_var then _sel_n = 0;

  _keep = 1;
  do _i = 1 to _sel_n;
    if (&pos_var >= _sel_dis_st[_i]) and (&pos_var <= _sel_dis_end[_i]) then do;
      _keep = 0;
      leave;
    end;
  end;

  if _keep then do;
    output;
    _sel_n + 1;
    if _sel_n > dim(_sel_dis_st) then do;
      put "ERROR: get_top_signal_within_dist exceeded the temporary selection buffer size of 50000.";
      stop;
    end;
    _sel_dis_st[_sel_n] = dis_st;
    _sel_dis_end[_sel_n] = dis_end;
  end;

  drop _sel_n _keep _i;
run;

%if &select_smallest_signal=1 %then %do;
proc sort data=&dsdout;
  by &grp_var &signal_var &pos_var;
run;
%end;
%else %do;
proc sort data=&dsdout;
  by &grp_var descending &signal_var &pos_var;
run;
%end;

proc datasets library=work nolist;
  delete _top_signal_candidates;
quit;

%mend;

/*Demo;

data tops;
input chr $ P BP;
cards;
2 0.001 1
2 0.011 40
2 0.00001 1000
2 0.0000001 1002
2 0.0000001 1002
2 0.0005 1500
2 0.001 3000000
2 0.001 90
2 0.01 40000
2 0.00001 10000
2 0.0005 150000
2 0.001 300000000
;
run;

options mprint mlogic symbolgen;

%get_top_signal_within_dist(dsdin=tops
                           ,grp_var=chr
                           ,signal_var=P
                           ,select_smallest_signal=1
                           ,pos_var=BP
                           ,pos_dist_thrshd=1000000
                           ,dsdout=tops1
                           ,signal_thrshd=1e-3);
*/
