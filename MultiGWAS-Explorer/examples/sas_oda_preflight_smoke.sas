options notes;
%macro quoted_token;
data preflight_token;
    length token $20;
    token='/*';
run;
/* An actual comment follows the quoted delimiter. */
%mend;
%quoted_token;

data preflight_value; * author's explanatory note; value=62; run;
/* Documentation example that must not be treated as executed SAS:
%put training; explanatory text;
*/
data _null_;
    set preflight_token;
    if token='/*' then put 'PREFLIGHT_TOKEN_OK';
    else put 'ERROR: preflight token was not preserved';
run;
data _null_;
    set preflight_value;
    if value=62 then put 'PREFLIGHT_VALUE_OK';
    else put 'ERROR: inline statement comment was not accepted';
run;
proc print data=preflight_token; run;
proc print data=preflight_value; run;
