/* Output-copy helper; it does not recalculate test statistics. */
%macro af_export(data=,name=);
  data work.af_export_copy;set &data;format _numeric_ best32.;run;
  proc export data=work.af_export_copy
    outfile="&AF_OUT/&name..csv" dbms=csv replace;
  run;
%mend;
