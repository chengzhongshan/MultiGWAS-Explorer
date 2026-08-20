%include "~/QueryLD_SNPs_at_Haploreg4.sas";
%include "~/QueryMulti_LD_SNPs_at_Haploreg4.sas";

%QueryMulti_LD_SNPs_at_Haploreg4(
  snps=rs185665940 rs10166057,
  LDpop=EUR,
  ldThresh=0.5,
  dsd4LD_SNPs=LD_SNPs_EUR
);

%QueryMulti_LD_SNPs_at_Haploreg4(
  snps=rs185665940 rs10166057,
  LDpop=ASN,
  ldThresh=0.5,
  dsd4LD_SNPs=LD_SNPs_ASN
);

proc export data=LD_SNPs_EUR outfile="~/reviewer_haploreg_ld_eur_r2ge0p5.tsv"
  dbms=tab replace;
run;
proc export data=LD_SNPs_ASN outfile="~/reviewer_haploreg_ld_asn_r2ge0p5.tsv"
  dbms=tab replace;
run;

proc contents data=LD_SNPs_EUR; run;
proc print data=LD_SNPs_EUR(obs=20); run;
proc contents data=LD_SNPs_ASN; run;
proc print data=LD_SNPs_ASN(obs=20); run;
