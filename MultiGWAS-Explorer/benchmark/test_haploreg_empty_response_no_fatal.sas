/* Regression: an expected empty HaploReg response must not emit SAS ERROR. */
%include "~/QueryLD_SNPs_at_Haploreg4.sas";

%QueryLD_SNPs_at_Haploreg4(
  snp=rs113595177,
  LDpop=EUR,
  ldThresh=0.2,
  outdsd=empty_response_test
);

data _null_;
  put "NOTE: EMPTY_RESPONSE_REGRESSION query_ok=&HAPLOREG_LD_QUERY_OK"
      " source=&HAPLOREG_LD_QUERY_SOURCE http=&HAPLOREG_LD_HTTP_STATUS";
run;
