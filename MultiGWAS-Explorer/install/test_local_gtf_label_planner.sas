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

data work.gtf_label_cd55_smoke;
  length SNP $20;
  input SNP $ BP;
datalines;
rs2802216 207394658
rs200078770 207394677
rs4303075 207394678
rs11117664 207394689
;
run;

%Plan_Local_GTF_Target_Labels(
  wide_dsd=work.gtf_label_cd55_smoke,
  label_snps=rs2802216 rs200078770 rs4303075 rs11117664,
  center_bp=207394658,
  window_bp=1000031,
  design_width=950,
  design_height=1000,
  font_size=10,
  out_angle=smoke_angle,
  out_positions=smoke_positions,
  out_headroom=smoke_headroom,
  out_center_offset=smoke_offset
);
%if &smoke_angle ne 0 or %sysfunc(countw(%superq(smoke_positions),%str( ))) ne 4
    or %sysevalf(&smoke_headroom>=0.05) %then %do;
  %put ERROR: SAS four-target local-GTF label planner smoke test failed.;
  %abort 255;
%end;

data work.gtf_label_chr21_smoke;
  length SNP $20;
  input SNP $ BP;
datalines;
rs2070788 41470061
rs383510 41486440
rs462687 41425130
rs429442 41489405
;
run;

%Plan_Local_GTF_Target_Labels(
  wide_dsd=work.gtf_label_chr21_smoke,
  label_snps=rs2070788 rs383510 rs462687 rs429442,
  center_bp=41470061,
  window_bp=1044931,
  design_width=950,
  design_height=1000,
  font_size=10,
  out_angle=smoke_angle,
  out_positions=smoke_positions,
  out_headroom=smoke_headroom,
  out_center_offset=smoke_offset
);
%if &smoke_angle ne 0 or %sysfunc(countw(%superq(smoke_positions),%str( ))) ne 4
    or %sysevalf(&smoke_headroom>=0.05) %then %do;
  %put ERROR: SAS four-target chr21 horizontal-label smoke test failed.;
  %abort 255;
%end;

data work.gtf_label_crowded_smoke;
  length SNP $20;
  input SNP $ BP;
datalines;
rs1234567891 207394650
rs1234567892 207394651
rs1234567893 207394652
rs1234567894 207394653
rs1234567895 207394654
rs1234567896 207394655
rs1234567897 207394656
rs1234567898 207394657
;
run;

%Plan_Local_GTF_Target_Labels(
  wide_dsd=work.gtf_label_crowded_smoke,
  label_snps=rs1234567891 rs1234567892 rs1234567893 rs1234567894
             rs1234567895 rs1234567896 rs1234567897 rs1234567898,
  center_bp=207394650,
  window_bp=1000000,
  design_width=950,
  design_height=1000,
  font_size=10,
  out_angle=smoke_angle,
  out_positions=smoke_positions,
  out_headroom=smoke_headroom,
  out_center_offset=smoke_offset
);
%if &smoke_angle ne 90 or %sysfunc(countw(%superq(smoke_positions),%str( ))) ne 8
    or %sysevalf(&smoke_headroom<=0.1) %then %do;
  %put ERROR: SAS width-limited vertical-label smoke test failed.;
  %abort 255;
%end;
