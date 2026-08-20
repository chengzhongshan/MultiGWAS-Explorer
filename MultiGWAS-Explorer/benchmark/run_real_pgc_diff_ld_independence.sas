/* LD-clump the MAF-filtered real PGC schizophrenia differential candidates. */
options validvarname=any;

%include "~/QueryLD_SNPs_at_Haploreg4.sas";
%include "~/get_top_signal_with_ld.sas";

proc import
  datafile="~/PGC_SCZ_SAS_local_top_hits_manhattan_common_top_hits.csv"
  out=pgc_diff_candidates
  dbms=csv replace;
  guessingrows=max;
run;

%get_top_signal_with_ld(
  dsdin=pgc_diff_candidates,
  snp_var=SNP,
  grp_var=CHR,
  signal_var=ALL_DIFF_P,
  select_smallest_signal=1,
  pos_var=BP,
  signal_thrshd=1e-6,
  ld_pops=EUR ASN,
  ld_r2_threshold=0.1,
  ld_population_rule=ANY,
  query_failure_action=DISTANCE,
  fallback_dist_bp=1e6,
  dsdout=pgc_diff_ld_leads,
  audit_out=pgc_diff_ld_audit
);

proc export data=pgc_diff_ld_leads
  outfile="~/reviewer_pgc_diff_ld_independent_leads.tsv"
  dbms=tab replace;
run;
proc export data=pgc_diff_ld_audit
  outfile="~/reviewer_pgc_diff_ld_audit.tsv"
  dbms=tab replace;
run;

proc print data=pgc_diff_ld_leads; run;
proc print data=pgc_diff_ld_audit; run;
