/* Regression test for the SAS ODA low-WORK-space local-GTF update path. */
options notes;

data scz_mh;
  input CHR BP SNP $ P1 P2;
  datalines;
21 100 rs1 0.01 0.20
21 200 rs2 0.03 0.04
;
run;

proc datasets library=work nolist;
  change scz_mh=_scz_mh_common_base;
quit;

data scz_mh / view=scz_mh;
  set _scz_mh_common_base;
  COMMON_ASSOC_P=min(of P1 P2);
run;

data final;
  length target $25;
  input target $ pos newpos top_y4label subgroup y1 _pos1_ pfactor;
  datalines;
rs1 100 100 9 1 -1 100 2
rs1 101 101 9 1 -1 100 2
rs2 200 200 9 1 -2 200 2
.   300 300 9 1 -3 300 2
;
run;

data final;
  if _n_=1 then do;
    declare hash seen();
    seen.defineKey('target');
    seen.defineData('target');
    seen.defineDone();
  end;
  modify final;
  if target ne '' and pos ne . then do;
    if seen.check()=0 then target='';
    else seen.add();
  end;
  replace;
run;

data _xtag_;
  length target $25;
  input target $ pos newpos;
  datalines;
rs1 100 110
rs2 200 190
;
run;

data final;
  if _n_=1 then do;
    declare hash label_positions(dataset:'_xtag_', duplicate:'r');
    label_positions.defineKey('target', 'pos');
    label_positions.defineData('newpos');
    label_positions.defineDone();
  end;
  modify final;
  if target ne '' and pos ne . then do;
    if label_positions.find()=0 then subgroup=0;
  end;
  replace;
run;

data final_gene_lookup;
  input y1 _pos1_ pfactor;
  datalines;
-1 100 90
-2 200 210
;
run;

data final;
  if _n_=1 then do;
    declare hash gene_positions(dataset:'final_gene_lookup', duplicate:'r');
    gene_positions.defineKey('y1', '_pos1_');
    gene_positions.defineData('pfactor');
    gene_positions.defineDone();
  end;
  modify final;
  pfactor=.;
  if gene_positions.find()=0 and pfactor ne . then _pos1_=pfactor;
  replace;
run;

proc sql noprint;
  select count(*) into :label_count trimmed from final where target ne '';
  select count(*) into :adjusted_label_count trimmed
    from final where (target='rs1' and newpos=110) or (target='rs2' and newpos=190);
  select count(*) into :adjusted_gene_count trimmed
    from final where (y1=-1 and _pos1_=90) or (y1=-2 and _pos1_=210);
  select count(*) into :common_view_count trimmed
    from scz_mh where COMMON_ASSOC_P in (0.01,0.03);
quit;

%macro assert_equal(actual, expected, name);
  %if &actual ne &expected %then %do;
    %put ERROR: LOW_WORK_TEST &name expected=&expected actual=&actual;
    %abort cancel;
  %end;
%mend;

%assert_equal(&label_count,2,label_count);
%assert_equal(&adjusted_label_count,2,adjusted_label_count);
%assert_equal(&adjusted_gene_count,3,adjusted_gene_count);
%assert_equal(&common_view_count,2,common_view_count);
%put NOTE: LOW_WORK_HASH_UPDATE_TEST_PASSED;
