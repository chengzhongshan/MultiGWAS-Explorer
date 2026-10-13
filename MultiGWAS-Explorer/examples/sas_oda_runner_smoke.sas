/* Small file-based smoke test; expected n=3 and sum=60. */
data runner_smoke;
    input id value;
    datalines;
1 10
2 20
3 30
;
run;
title "SAS ODA runner smoke test";
proc means data=runner_smoke n sum;
    var value;
run;
proc print data=runner_smoke; run;
title;
