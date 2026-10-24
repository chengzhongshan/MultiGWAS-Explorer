%macro rows_in_sas_dsd(test_dsd);
%local dsid nlobs rc;
%let nlobs=0;
%let dsid=%sysfunc(open(&test_dsd,IS));
%if &dsid=0 %then 
%put %sysfunc(sysmsg());

%let nlobs=%sysfunc(attrn(&dsid,NLOBS));
%let rc=%sysfunc(close(&dsid));

&nlobs

%if &nlobs gt 1 %then %do;
 %put There are &nlobs in your dataset;
%end;

%mend;
