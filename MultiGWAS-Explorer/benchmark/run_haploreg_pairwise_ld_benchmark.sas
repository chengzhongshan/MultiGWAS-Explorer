%include "~/QueryLD_SNPs_at_Haploreg4.sas";

%QueryLD_SNPs_at_Haploreg4(
  snp=rs185665940, LDpop=EUR, ldThresh=0.01, outdsd=LD_PAIR_EUR_185
);
%QueryLD_SNPs_at_Haploreg4(
  snp=rs185665940, LDpop=ASN, ldThresh=0.01, outdsd=LD_PAIR_ASN_185
);
%QueryLD_SNPs_at_Haploreg4(
  snp=rs10166057, LDpop=EUR, ldThresh=0.01, outdsd=LD_PAIR_EUR_101
);
%QueryLD_SNPs_at_Haploreg4(
  snp=rs10166057, LDpop=ASN, ldThresh=0.01, outdsd=LD_PAIR_ASN_101
);

data LD_PAIR_EUR_185_C; set LD_PAIR_EUR_185; keep chr--EUR query_snp_rsid; run;
data LD_PAIR_ASN_185_C; set LD_PAIR_ASN_185; keep chr--EUR query_snp_rsid; run;
data LD_PAIR_EUR_101_C; set LD_PAIR_EUR_101; keep chr--EUR query_snp_rsid; run;
data LD_PAIR_ASN_101_C; set LD_PAIR_ASN_101; keep chr--EUR query_snp_rsid; run;

data PAIRWISE_LD;
  length population $3 source_query $20;
  set LD_PAIR_EUR_185_C(in=a) LD_PAIR_ASN_185_C(in=b)
      LD_PAIR_EUR_101_C(in=c) LD_PAIR_ASN_101_C(in=d);
  if a then do; population='EUR'; source_query='rs185665940'; end;
  else if b then do; population='ASN'; source_query='rs185665940'; end;
  else if c then do; population='EUR'; source_query='rs10166057'; end;
  else if d then do; population='ASN'; source_query='rs10166057'; end;
  if rsID in ('rs185665940','rs10166057');
  keep population source_query rsID r2 D: is_query_snp chr pos_hg38 ref alt AFR AMR ASN EUR;
run;

proc export data=PAIRWISE_LD
  outfile="~/haploreg_pairwise_ld.tsv"
  dbms=tab replace;
run;
proc print data=PAIRWISE_LD; run;
