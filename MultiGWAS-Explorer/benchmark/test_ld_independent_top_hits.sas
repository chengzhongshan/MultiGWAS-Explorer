/* Reviewer-facing regression fixture: LD clumping must retain two independent
   signals that a 1 Mb physical-distance rule collapses into one locus. */
options validvarname=any;

%include "~/QueryLD_SNPs_at_Haploreg4.sas";
%include "~/get_top_signal_with_ld.sas";
%include "~/get_top_signal_within_dist.sas";

data fixture_candidates;
  length SNP $20;
  input CHR BP SNP $ P;
cards;
2 72142747 rs10166057  1e-9
2 72134376 rs3768644   2e-9
2 72041898 rs185665940 3e-9
2 72023810 rs139205793 4e-9
;
run;

%get_top_signal_with_ld(
  dsdin=fixture_candidates,
  snp_var=SNP,
  grp_var=CHR,
  signal_var=P,
  select_smallest_signal=1,
  pos_var=BP,
  signal_thrshd=1e-6,
  ld_pops=EUR,
  ld_r2_threshold=0.5,
  ld_population_rule=ANY,
  query_failure_action=DISTANCE,
  fallback_dist_bp=1e6,
  dsdout=fixture_ld_leads,
  audit_out=fixture_ld_audit
);

%get_top_signal_within_dist(
  dsdin=fixture_candidates,
  grp_var=CHR,
  signal_var=P,
  select_smallest_signal=1,
  pos_var=BP,
  pos_dist_thrshd=1e6,
  dsdout=fixture_distance_leads,
  signal_thrshd=1e-6
);

proc export data=fixture_ld_leads
  outfile="~/reviewer_fixture_ld_independent_leads.tsv"
  dbms=tab replace;
run;
proc export data=fixture_ld_audit
  outfile="~/reviewer_fixture_ld_audit.tsv"
  dbms=tab replace;
run;
proc export data=fixture_distance_leads
  outfile="~/reviewer_fixture_distance_leads.tsv"
  dbms=tab replace;
run;

proc print data=fixture_ld_leads; run;
proc print data=fixture_ld_audit; run;
proc print data=fixture_distance_leads; run;
