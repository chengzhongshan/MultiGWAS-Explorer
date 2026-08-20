/*
Greedy, P-value-ranked LD clumping for common or differential GWAS hits.

The best remaining candidate is selected as a lead. HaploReg is queried for
that lead in each requested ancestry, and candidate SNPs with r2 at or above
ld_r2_threshold are removed. Physical distance is used only when every LD
query for a lead fails (or the lead is not an rsID), and every selected/pruned
candidate is recorded in audit_out.
*/
%macro get_top_signal_with_ld(
  dsdin=,
  snp_var=SNP,
  grp_var=CHR,
  signal_var=P,
  select_smallest_signal=1,
  pos_var=BP,
  signal_thrshd=1,
  ld_pops=EUR,
  ld_r2_threshold=0.1,
  ld_population_rule=ANY,
  query_failure_action=DISTANCE,
  fallback_dist_bp=1e6,
  max_leads=0,
  dsdout=top_ind,
  audit_out=top_hit_ld_audit
);
%local _candidate_n _iter _remaining_n _lead_id _lead_snp _lead_chr
       _lead_query_snp _lead_bp _lead_signal _lead_is_rs _pop_n _pi _pop _valid_pop_n
       _requested_pop_n _response_valid _response_self_n _selected_n
       _rule _failure_action _max_leads_num _has_rsid _has_r2 _cache_pop_n
       _lead_query_status;
%global TOP_HIT_LD_SELECTION_STATUS TOP_HIT_LD_SELECTED_N TOP_HIT_LD_PRUNED_N;

%let TOP_HIT_LD_SELECTION_STATUS=RUNNING;
%let TOP_HIT_LD_SELECTED_N=0;
%let TOP_HIT_LD_PRUNED_N=0;
%let _rule=%upcase(%sysfunc(strip(%superq(ld_population_rule))));
%let _failure_action=%upcase(%sysfunc(strip(%superq(query_failure_action))));
%let _requested_pop_n=%sysfunc(countw(%superq(ld_pops),%str( )));
%let _max_leads_num=%sysfunc(inputn(%superq(max_leads),best32.));
%if %sysevalf(%superq(_max_leads_num)=,boolean) %then %let _max_leads_num=0;

%if &_requested_pop_n=0 %then %do;
  %put ERROR: get_top_signal_with_ld requires at least one HaploReg population in ld_pops=.;
  %let TOP_HIT_LD_SELECTION_STATUS=INVALID_LD_POPULATIONS;
  %return;
%end;
%if &_rule ne ANY and &_rule ne ALL %then %do;
  %put ERROR: ld_population_rule must be ANY or ALL, not &ld_population_rule.;
  %let TOP_HIT_LD_SELECTION_STATUS=INVALID_POPULATION_RULE;
  %return;
%end;

data _ld_candidates_unsorted;
  set &dsdin;
  length _ld_snp_key $128 _ld_chr_key $32;
  length _ld_signal _ld_bp_num 8;
  _ld_snp_key=upcase(strip(vvalue(&snp_var)));
  _ld_chr_key=upcase(compress(strip(vvalue(&grp_var)),' '));
  _ld_chr_key=prxchange('s/^CHR//i',1,_ld_chr_key);
  _ld_bp_num=input(strip(vvalue(&pos_var)),best32.);
  _ld_signal=input(strip(vvalue(&signal_var)),best32.);
  if missing(_ld_snp_key) or missing(_ld_chr_key) or missing(_ld_bp_num) then delete;
  if missing(_ld_signal) or _ld_signal <= 0 or _ld_signal > &signal_thrshd then delete;
run;

proc sort data=_ld_candidates_unsorted out=_ld_candidates_ranked;
%if &select_smallest_signal=1 %then %do;
  by _ld_signal _ld_chr_key _ld_bp_num _ld_snp_key;
%end;
%else %do;
  by descending _ld_signal _ld_chr_key _ld_bp_num _ld_snp_key;
%end;
run;

data _ld_candidates;
  if _n_=1 then do;
    declare hash _seen();
    _seen.defineKey('_ld_snp_key');
    _seen.defineDone();
  end;
  set _ld_candidates_ranked;
  if _seen.check()=0 then delete;
  _seen.add();
  length _ld_pruned_by $128 _ld_prune_population $32
         _ld_query_status $24 _ld_independence_method $24;
  retain _ld_candidate_id 0;
  _ld_candidate_id+1;
  _ld_selected=0;
  _ld_pruned=0;
  _ld_lead_rank=.;
  _ld_pruned_lead_rank=.;
  _ld_prune_r2=.;
run;

proc sql noprint;
  select count(*) into :_candidate_n trimmed from _ld_candidates;
quit;
%if %sysevalf(%superq(_candidate_n)=,boolean) %then %let _candidate_n=0;

%if &_candidate_n=0 %then %do;
  data &dsdout; set _ld_candidates; run;
  data &audit_out;
    length candidate_snp lead_snp $128 selection_action $24;
    stop;
  run;
  %let TOP_HIT_LD_SELECTION_STATUS=NO_CANDIDATES;
  %goto _ld_cleanup;
%end;

%do _iter=1 %to &_candidate_n;
  %let _remaining_n=0;
  %let _lead_id=;
  %let _lead_snp=;
  %let _lead_query_snp=;
  %let _lead_chr=;
  %let _lead_bp=;
  %let _lead_signal=;
  proc sql noprint;
    select count(*) into :_remaining_n trimmed
    from _ld_candidates
    where _ld_selected=0 and _ld_pruned=0;

    select _ld_candidate_id, _ld_snp_key, _ld_chr_key, _ld_bp_num, _ld_signal
      into :_lead_id trimmed, :_lead_snp trimmed, :_lead_chr trimmed,
           :_lead_bp trimmed, :_lead_signal trimmed
    from _ld_candidates
    where _ld_selected=0 and _ld_pruned=0
    having _ld_candidate_id=min(_ld_candidate_id);
  quit;

  %if %sysevalf(%superq(_remaining_n)=,boolean) %then %let _remaining_n=0;
  %if &_remaining_n=0 or %sysevalf(%superq(_lead_id)=,boolean) %then %goto _ld_finished;

  proc sql noprint;
    select count(*) into :_selected_n trimmed
    from _ld_candidates where _ld_selected=1;
  quit;
  %if %sysevalf(%superq(_selected_n)=,boolean) %then %let _selected_n=0;
  %if %sysevalf(&_max_leads_num > 0) and %sysevalf(&_selected_n >= &_max_leads_num) %then %goto _ld_finished;

  proc sql;
    update _ld_candidates
    set _ld_selected=1,
        _ld_lead_rank=&_selected_n+1
    where _ld_candidate_id=&_lead_id;
  quit;

  data _null_;
    call symputx('_lead_is_rs',prxmatch('/^RS[0-9]+$/i',symget('_lead_snp'))>0,'L');
  run;

  data _ld_lead_edges;
    length query_snp proxy_snp $128 ld_population $16 proxy_r2 8;
    stop;
  run;
  %let _valid_pop_n=0;
  %let _cache_pop_n=0;

  %if &_lead_is_rs %then %do;
    %let _lead_query_snp=%sysfunc(lowcase(%superq(_lead_snp)));
    %let _pop_n=%sysfunc(countw(%superq(ld_pops),%str( )));
    %do _pi=1 %to &_pop_n;
      %let _pop=%scan(%superq(ld_pops),&_pi,%str( ));
      %QueryLD_SNPs_at_Haploreg4(
        snp=&_lead_query_snp,
        LDpop=&_pop,
        ldThresh=&ld_r2_threshold,
        outdsd=_ld_raw_response
      );

      %let _has_rsid=0;
      %let _has_r2=0;
      %if %sysfunc(exist(work._ld_raw_response)) %then %do;
        proc sql noprint;
          select sum(upcase(name)='RSID'), sum(upcase(name)='R2')
            into :_has_rsid trimmed, :_has_r2 trimmed
          from dictionary.columns
          where libname='WORK' and memname='_LD_RAW_RESPONSE';
        quit;
      %end;
      %if %sysevalf(%superq(_has_rsid)=,boolean) %then %let _has_rsid=0;
      %if %sysevalf(%superq(_has_r2)=,boolean) %then %let _has_r2=0;

      %if &_has_rsid>0 and &_has_r2>0 %then %do;
        data _ld_response_norm;
          length query_snp proxy_snp $128 ld_population $16 proxy_r2 8;
          set _ld_raw_response;
          query_snp="&_lead_snp";
          proxy_snp=upcase(strip(vvaluex('rsID')));
          ld_population="%upcase(&_pop)";
          proxy_r2=input(strip(vvaluex('r2')),best32.);
          if not missing(proxy_snp) and not missing(proxy_r2);
          keep query_snp proxy_snp ld_population proxy_r2;
        run;
        proc sql noprint;
          select count(*) into :_response_self_n trimmed
          from _ld_response_norm
          where proxy_snp="&_lead_snp";
        quit;
        %if %sysevalf(%superq(_response_self_n)=,boolean) %then %let _response_self_n=0;
        %if &_response_self_n>0 %then %do;
          %let _valid_pop_n=%eval(&_valid_pop_n+1);
          %if %upcase(%superq(HAPLOREG_LD_QUERY_SOURCE))=LOCAL_CACHE %then
            %let _cache_pop_n=%eval(&_cache_pop_n+1);
          proc append base=_ld_lead_edges data=_ld_response_norm force; run;
        %end;
      %end;
      proc datasets library=work nolist;
        delete _ld_raw_response _ld_response_norm;
      quit;
    %end;
  %end;

  %if &_valid_pop_n>0 %then %do;
    %if &_valid_pop_n=&_requested_pop_n %then %let _lead_query_status=OK;
    %else %let _lead_query_status=PARTIAL;
    %if &_cache_pop_n=&_valid_pop_n %then %let _lead_query_status=&_lead_query_status._LOCAL_CACHE;
    %else %if &_cache_pop_n>0 %then %let _lead_query_status=&_lead_query_status._MIXED;

    proc sql;
      create table _ld_prune_keys as
      select proxy_snp as _ld_snp_key length=128,
             max(proxy_r2) as _edge_r2,
             case
               when count(distinct ld_population)>1 then 'MULTI'
               else min(ld_population)
             end as _edge_population length=32
      from _ld_lead_edges
      where proxy_snp ne "&_lead_snp"
        and proxy_r2 >= &ld_r2_threshold
      group by proxy_snp
      %if &_rule=ALL %then %do;
      having count(distinct ld_population) >= &_requested_pop_n
      %end;
      ;
    quit;

    proc sort data=_ld_candidates; by _ld_snp_key; run;
    proc sort data=_ld_prune_keys; by _ld_snp_key; run;
    data _ld_candidates;
      merge _ld_candidates(in=_candidate) _ld_prune_keys;
      by _ld_snp_key;
      if _candidate;
      if _ld_selected=0 and _ld_pruned=0 and not missing(_edge_r2) then do;
        _ld_pruned=1;
        _ld_pruned_by="&_lead_snp";
        _ld_pruned_lead_rank=&_selected_n+1;
        _ld_prune_r2=_edge_r2;
        _ld_prune_population=_edge_population;
        _ld_independence_method='LD';
      end;
      drop _edge_r2 _edge_population;
    run;
    proc sort data=_ld_candidates; by _ld_candidate_id; run;

    proc sql;
      update _ld_candidates
      set _ld_query_status="&_lead_query_status",
          _ld_independence_method='LD'
      where _ld_candidate_id=&_lead_id;
    quit;
  %end;
  %else %do;
    %if &_failure_action=DISTANCE and %sysevalf(&fallback_dist_bp > 0) %then %do;
      data _ld_candidates;
        set _ld_candidates;
        if _ld_selected=0 and _ld_pruned=0
           and _ld_chr_key="&_lead_chr"
           and _ld_bp_num >= (&_lead_bp-0.5*&fallback_dist_bp)
           and _ld_bp_num <= (&_lead_bp+0.5*&fallback_dist_bp) then do;
          _ld_pruned=1;
          _ld_pruned_by="&_lead_snp";
          _ld_pruned_lead_rank=&_selected_n+1;
          _ld_prune_population='NA';
          _ld_independence_method='DISTANCE_FALLBACK';
        end;
        if _ld_candidate_id=&_lead_id then do;
          _ld_query_status='NO_LD_RESPONSE';
          _ld_independence_method='DISTANCE_FALLBACK';
        end;
      run;
    %end;
    %else %do;
      proc sql;
        update _ld_candidates
        set _ld_query_status='NO_LD_RESPONSE',
            _ld_independence_method='UNRESOLVED_NO_LD'
        where _ld_candidate_id=&_lead_id;
      quit;
    %end;
  %end;

  proc datasets library=work nolist;
    delete _ld_lead_edges _ld_prune_keys;
  quit;
%end;

%_ld_finished:
data &dsdout;
  set _ld_candidates;
  where _ld_selected=1;
  length INDEPENDENCE_METHOD $24 LD_POPULATIONS $64 LD_POPULATION_RULE $8
         LD_QUERY_STATUS $24;
  LD_LEAD_RANK=_ld_lead_rank;
  INDEPENDENCE_METHOD=_ld_independence_method;
  LD_R2_THRESHOLD=&ld_r2_threshold;
  LD_POPULATIONS="%upcase(&ld_pops)";
  LD_POPULATION_RULE="&_rule";
  LD_QUERY_STATUS=_ld_query_status;
  LD_FALLBACK_BP=&fallback_dist_bp;
  drop _ld_:;
run;
proc sort data=&dsdout; by LD_LEAD_RANK; run;

data &audit_out;
  set _ld_candidates;
  length candidate_snp lead_snp $128 selection_action $24
         independence_method $24 query_status $24 ld_population $32;
  candidate_rank=_ld_candidate_id;
  candidate_snp=_ld_snp_key;
  chr=_ld_chr_key;
  bp=_ld_bp_num;
  signal=_ld_signal;
  if _ld_selected then do;
    selection_action='SELECTED_LEAD';
    lead_rank=_ld_lead_rank;
    lead_snp=_ld_snp_key;
  end;
  else if _ld_pruned then do;
    selection_action=ifc(_ld_independence_method='LD','PRUNED_LD','PRUNED_DISTANCE_FALLBACK');
    lead_rank=_ld_pruned_lead_rank;
    lead_snp=_ld_pruned_by;
  end;
  else selection_action='NOT_REACHED';
  prune_r2=_ld_prune_r2;
  ld_population=_ld_prune_population;
  query_status=_ld_query_status;
  independence_method=_ld_independence_method;
  keep candidate_rank candidate_snp chr bp signal selection_action lead_rank
       lead_snp prune_r2 ld_population query_status independence_method;
run;
proc sort data=&audit_out; by candidate_rank; run;

proc sql noprint;
  select count(*) into :TOP_HIT_LD_SELECTED_N trimmed
  from _ld_candidates where _ld_selected=1;
  select count(*) into :TOP_HIT_LD_PRUNED_N trimmed
  from _ld_candidates where _ld_pruned=1;
quit;
%let TOP_HIT_LD_SELECTION_STATUS=COMPLETE;
%put NOTE: LD top-hit selection complete: candidates=&_candidate_n selected=&TOP_HIT_LD_SELECTED_N pruned=&TOP_HIT_LD_PRUNED_N r2=&ld_r2_threshold populations=&ld_pops rule=&_rule.;

%_ld_cleanup:
proc datasets library=work nolist;
  delete _ld_candidates_unsorted _ld_candidates_ranked _ld_candidates
         _ld_lead_edges _ld_prune_keys _ld_raw_response _ld_response_norm;
quit;
%mend;
