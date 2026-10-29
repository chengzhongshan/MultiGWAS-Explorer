%macro sas_dsd_exist(sasdsd,lib=work);
%local sas_dt_exist;
%let sas_dt_exist=0;
%if %sysfunc(prxmatch(/\./,&sasdsd)) %then %do;
   %let lib=%scan(&sasdsd,1,.);
   %let sasdsd=%scan(&sasdsd,2,.);
%end;
%let sas_dt_exist=%sysfunc(exist(&lib..&sasdsd));
&sas_dt_exist
%mend;
