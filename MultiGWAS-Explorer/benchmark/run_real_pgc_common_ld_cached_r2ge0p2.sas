/* Fast common-association clumping from the official HaploReg v4 LD archives. */
options validvarname=any;

%include "~/QueryLD_SNPs_at_Haploreg4.sas";
%include "~/get_top_signal_with_ld.sas";

proc import
  datafile="~/pgc_common_ld_candidates.csv"
  out=pgc_common_candidates
  dbms=csv replace;
  guessingrows=max;
run;

proc import
  datafile="~/pgc_common_haploreg_candidate_ld.tsv"
  out=pgc_common_haploreg_cache
  dbms=tab replace;
  guessingrows=max;
run;
proc datasets library=work nolist;
  modify pgc_common_haploreg_cache;
  index create query_population=(query_snp ld_population);
quit;

%let HAPLOREG_LD_CACHE_DSD=work.pgc_common_haploreg_cache;
%let HAPLOREG_LD_CACHE_MIN_R2=0.2;
%let HAPLOREG_LD_CACHE_MISS_ACTION=WEB;

data pgc_common_candidates;
  set pgc_common_candidates;
  COMMON_ASSOC_P=min(
    ALL_GROUP1_P, ALL_GROUP2_P,
    EUR_GROUP1_P, EUR_GROUP2_P,
    ASN_GROUP1_P, ASN_GROUP2_P
  );
run;

%get_top_signal_with_ld(
  dsdin=pgc_common_candidates,
  snp_var=SNP,
  grp_var=CHR,
  signal_var=COMMON_ASSOC_P,
  select_smallest_signal=1,
  pos_var=BP,
  signal_thrshd=5e-8,
  ld_pops=EUR ASN,
  ld_r2_threshold=0.2,
  ld_population_rule=ANY,
  query_failure_action=DISTANCE,
  fallback_dist_bp=1e6,
  dsdout=pgc_common_ld_leads,
  audit_out=pgc_common_ld_audit
);

proc export data=pgc_common_ld_leads
  outfile="~/pgc_common_ld_cached_r2ge0p2_leads.tsv"
  dbms=tab replace;
run;
proc export data=pgc_common_ld_audit
  outfile="~/pgc_common_ld_cached_r2ge0p2_audit.tsv"
  dbms=tab replace;
run;

proc freq data=pgc_common_ld_audit;
  tables selection_action independence_method query_status / missing;
run;
proc print data=pgc_common_ld_leads(obs=25); run;
