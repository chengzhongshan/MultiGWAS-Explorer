%macro check_var_type(dsdin,var_rgx);
proc contents data=&dsdin out=&dsdin._fmt noprint;
run;
proc sql noprint;
/*select distinct(TYPE) into: var_type separated by '|'*/
select TYPE into: var_type separated by '|'
from &dsdin._fmt
where prxmatch("/&var_rgx/i",NAME);

%put The types of variables matched with regular expression &var_rgx in dsd &dsdin are:;
%put &var_type;
%if %index(&var_type,|) %then %put Your regular expression matchs with >1 variables, please restrict your regular expression to only match with a unique variable name;

%mend;
