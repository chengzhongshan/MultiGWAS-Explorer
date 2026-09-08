
%macro ntokens(list);
    %if %length(%superq(list))=0 %then 0;
    %else %eval(1 + %length(%sysfunc(compbl(&list))) - %length(%sysfunc(compress(&list))));
%mend ntokens;
