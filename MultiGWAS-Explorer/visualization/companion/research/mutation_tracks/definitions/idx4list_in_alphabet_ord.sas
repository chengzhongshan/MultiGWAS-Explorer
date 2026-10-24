%macro idx4list_in_alphabet_ord(
/*This macro is very useful for obtaining numeric index in an alphabetic sorted list for a query list;
The output idx can be used in a macro loop or in a array at data step!*/
list= ,/*Blank space separated elements*/
outdsd=list_order, /*A sas dataset will be created to contain the index
for these elements sorted in alphabet mode*/
index_list_var=idx4list, /*a global macro var will be created to contain these index
Make sure the value for index_list_var is not index_list_var, as it will lead conflict between
the global macro var and the local macro var with the same name*/
sep=%nrbquote(,) /*separator for the newly created index*/
);
%global &index_list_var;
%let &index_list_var=;
%rank4grps(
grps=&list,
dsdout=&outdsd
);

proc sql noprint;
select char_ord into: &index_list_var separated by "&sep"
from &outdsd;

%mend;