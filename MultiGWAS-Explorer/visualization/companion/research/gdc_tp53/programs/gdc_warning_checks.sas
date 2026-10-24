/* Baseline warnings below are intentional evidence, not the corrected run. */
data work.fixture;
 length type $8 grp $8;
 input st end type $ grp $;
 datalines;
10 90 gene A
10 20 exon A
70 90 exon A
;
run;
data work.old_same;set work.fixture;run;
%BookBaselineGrouping(gene_bed_dsd=work.old_same,dsdout=work.old_same,gene_dist_thrhd=0);
%BookBaselineGrouping(gene_bed_dsd=work.fixture,dsdout=work.old_blank,reg_type=,gene_dist_thrhd=0);
%put GDC_FIXED_TESTS_BEGIN;
data work.fixed_same;set work.fixture;run;
%adj_grpnum4close_gene_bed_regs(gene_bed_dsd=work.fixed_same,dsdout=work.fixed_same,gene_dist_thrhd=0);
%adj_grpnum4close_gene_bed_regs(gene_bed_dsd=work.fixture,dsdout=work.fixed_blank,reg_type=,gene_dist_thrhd=0);
proc sort data=work.old_same;by _all_;run;
proc sort data=work.fixed_same;by _all_;run;
proc compare base=work.old_same compare=work.fixed_same noprint;run;
%let same_rc=&sysinfo;
proc sort data=work.old_blank;by _all_;run;
proc sort data=work.fixed_blank;by _all_;run;
proc compare base=work.old_blank compare=work.fixed_blank noprint;run;
%let blank_rc=&sysinfo;
data work.expected;length marker $8;marker='A';chr=17;pos=123;value=1;run;
data work.no_extra;set work.expected;run;
%BookGDCDropCoordinates(data=work.no_extra,chr=chr,pos=pos);
proc compare base=work.expected compare=work.no_extra noprint;run;
%let absent_rc=&sysinfo;
data work.has_extra;set work.expected;chr_BRCA=17;pos_BRCA=123;ChR_LUAD=17;POS_LUAD=123;run;
%BookGDCDropCoordinates(data=work.has_extra,chr=chr,pos=pos);
proc compare base=work.expected compare=work.has_extra noprint;run;
%let present_rc=&sysinfo;
data work.warning_checks;
 length case_name $32;
 case_name='group_same_input_output';compare_sysinfo=&same_rc;output;
 case_name='group_implicit_region_type';compare_sysinfo=&blank_rc;output;
 case_name='drop_absent_prefix';compare_sysinfo=&absent_rc;output;
 case_name='drop_existing_mixed_case_prefix';compare_sysinfo=&present_rc;output;
run;
%put GDC_FIXED_TESTS_END;
