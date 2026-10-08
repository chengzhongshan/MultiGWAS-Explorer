data work.gtf_label_planner_smoke;
  length SNP $20;
  input SNP $ BP;
datalines;
rs2070788 41400000
rs383510 41420000
;
run;

%include "~/Plan_Local_GTF_Target_Labels.sas";
%let smoke_angle=;
%let smoke_positions=;
%let smoke_headroom=;
%let smoke_offset=;
%Plan_Local_GTF_Target_Labels(
  wide_dsd=work.gtf_label_planner_smoke,
  label_snps=rs2070788 rs383510,
  center_bp=41400000,
  window_bp=1016379,
  design_width=950,
  design_height=1000,
  font_size=10,
  out_angle=smoke_angle,
  out_positions=smoke_positions,
  out_headroom=smoke_headroom,
  out_center_offset=smoke_offset
);
%put NOTE: LOCAL_GTF_LABEL_SMOKE angle=&smoke_angle positions=&smoke_positions headroom=&smoke_headroom offset=&smoke_offset;
%if &smoke_angle ne 0 or %length(%superq(smoke_positions))=0
    or %sysevalf(&smoke_headroom<=0) or %sysevalf(&smoke_offset<=0) %then %do;
  %put ERROR: SAS local-GTF label planner smoke test failed.;
  %abort 255;
%end;

%let horizontal_headroom=&smoke_headroom;
%Plan_Local_GTF_Target_Labels(
  wide_dsd=work.gtf_label_planner_smoke,
  label_snps=rs2070788 rs383510,
  center_bp=41400000,
  window_bp=1016379,
  design_width=950,
  design_height=1000,
  font_size=10,
  requested_layout=vertical,
  out_angle=smoke_angle,
  out_positions=smoke_positions,
  out_headroom=smoke_headroom,
  out_center_offset=smoke_offset
);
%if &smoke_angle ne 90 or %length(%superq(smoke_positions))=0
    or %sysevalf(&smoke_headroom<=&horizontal_headroom)
    or %sysevalf(&smoke_offset<=0) %then %do;
  %put ERROR: SAS vertical local-GTF label planner smoke test failed.;
  %abort 255;
%end;
