%macro QueryMulti_LD_SNPs_at_Haploreg4(
snps=rs17425819 rs2564978,/*Input multiple snps separated by blank space for querying at Haploreg4*/
LDpop=EUR,
ldThresh=0.5,
dsd4LD_SNPs=LD_SNPs /*Output a dataset into the work library;
Note: only common columns were kept, with SNP functional annotation
columns were removed due to different SNPs may have different
functional annotation columns and would fail in the union of these
datasets for all query SNPs; please use the macro QueryLD_SNPs_at_Haploreg4
to get full annotation table for each query SNP instead!*/
);
%local si LD_dsds dsd4LD_SNPs;
%let LD_dsds=;
%do si=1 %to %numargs(&snps);
%QueryLD_SNPs_at_Haploreg4(
snp=%scan(&snps,&si,%str( )),
LDpop=&LDpop,
ldThresh=&ldThresh,
outdsd=_LD_dsd&si
);
data _LD_dsd&si;
set _LD_dsd&si;
keep chr--EUR query_snp_rsid;
run;
%let LD_dsds=&LD_dsds _LD_dsd&si;
%end;

%Union_Data_In_Lib_Rgx(lib=work,
excluded=,
dsd_contain_rgx=_LD_dsd.*,
dsdout=&dsd4LD_SNPs
);
data &dsd4LD_SNPs;
set &dsd4LD_SNPs;
where rsID^="";
drop dsd;
run;
proc datasets nolist;
delete &LD_dsds;
run;
%mend;

/*Demo codes:;
%QueryLD_SNPs_at_Haploreg4(
snp=rs2564978,
LDpop=EUR,
ldThresh=0.95,
outdsd=LD_dsd
);

%QueryMulti_LD_SNPs_at_Haploreg4(
snps=rs17425819 rs2564978,
LDpop=EUR,
ldThresh=0.95,
dsd4LD_SNPs=LD_SNPs
);

*/
