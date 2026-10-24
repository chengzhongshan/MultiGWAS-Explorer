%MACRO GetVarAttrib(data, var, fun);
%local data_name dsid rc varnum attrib;

%let data_name = %scan(&data, 1, '(');

%* Initialize &attrib to empty in case there is an error;
%let attrib = ;
%* Open input dataset;
%let dsid = %sysfunc(open(&data_name));
	%** Note that data options are accepted in the open statement but I still use &data_name
	%** because I need it for the %PUT below in case there is an error;
%if &dsid = 0 %then
	%put GETVARATTRIB: ERROR - Dataset %upcase(&data_name) does not exist.;
%else %do;
	%let varnum = %sysfunc(varnum(&dsid,&var));
	%if &varnum = 0 %then
		%put GETVARATTRIB: ERROR - Variable %upcase(&var) does not exist in dataset %upcase(&data_name).;
	%else
		%let attrib = %sysfunc(&fun(&dsid,&varnum));
	%* Close input dataset;
	%let rc = %sysfunc(close(&dsid));
%end;

&attrib
%MEND GetVarAttrib;
