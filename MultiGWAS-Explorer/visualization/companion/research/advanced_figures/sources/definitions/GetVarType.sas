%MACRO GetVarType(data, var);
%local vartype;
%let vartype = %GetVarAttrib(&data, &var, vartype);
&vartype
%MEND GetVarType;
