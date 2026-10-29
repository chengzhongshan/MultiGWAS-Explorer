/* Drop only existing transpose-generated coordinate columns. */
%macro BookGDCDropCoordinates(data=,chr=,pos=);
 %local _id _n _i _name _drop _rc;
 %let _drop=;
 %let _id=%sysfunc(open(&data,i));
 %if &_id=0 %then %do;
  %put ERROR: BookGDCDropCoordinates cannot open &data.;
  %return;
 %end;
 %let _n=%sysfunc(attrn(&_id,nvars));
 %do _i=1 %to &_n;
  %let _name=%sysfunc(varname(&_id,&_i));
  %if %index(%upcase(%superq(_name)),%upcase(&chr._))=1 or
      %index(%upcase(%superq(_name)),%upcase(&pos._))=1 %then
   %let _drop=&_drop %sysfunc(nliteral(%superq(_name)));
 %end;
 %let _rc=%sysfunc(close(&_id));
 %if %length(%superq(_drop)) %then %do;
  data &data;set &data;drop &_drop;run;
 %end;
%mend;
