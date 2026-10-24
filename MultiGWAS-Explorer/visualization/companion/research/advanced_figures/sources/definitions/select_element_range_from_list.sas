%macro select_element_range_from_list(
list=,
st=1,	
end=,	/*If empty, the total number of elements will be used by it*/
sublist=sublist,/*Output a global macro variable for using by other SAS codes*/
sep=\s /*Perl regular expression for target separator (no | is allowed to match different separators)*/
);
%global &sublist;
%if %length(&end)=0 %then %let end=%ntokens(&list);
%let slcted_elems=%sysfunc(prxchange(s/^([^&sep]+&sep?){%sysevalf(&st-1)}(([^&sep]+&sep?){%sysevalf(&end-&st+1)}).*/\2/,-1,%bquote(&list))); 
*Need to remove the last separator;
*Note: the bquote for the macro var slcted_elems in case that it contain comma!;
%let slcted_elems=%sysfunc(prxchange(s/&sep$/ /,-1,%bquote(&slcted_elems)));
%put &slcted_elems;
%let &sublist=&slcted_elems;
%mend;
