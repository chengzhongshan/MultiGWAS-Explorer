/*
Repository-local HaploReg v4 LD query macro.

This is a hardened version of the project macro used by the reviewer
benchmark.  It preserves the original output shape and also publishes two
status macro variables so callers can distinguish an HTTP/import failure from
a successful LD response:

  HAPLOREG_LD_HTTP_STATUS
  HAPLOREG_LD_QUERY_OK
*/
%macro QueryLD_SNPs_at_Haploreg4(
  snp=rs2564978,
  LDpop=EUR,
  ldThresh=0.7,
  outdsd=LD_dsd
);
%local query_str _http_status _has_rsid _has_r2 _n_rows _out_lib _out_mem _response_bytes _response_records
       _cache_min _cache_action _cache_self_n;
%global HAPLOREG_LD_HTTP_STATUS HAPLOREG_LD_QUERY_OK HAPLOREG_LD_QUERY_SOURCE
        HAPLOREG_LD_CACHE_DSD HAPLOREG_LD_CACHE_MIN_R2 HAPLOREG_LD_CACHE_MISS_ACTION;
%let HAPLOREG_LD_HTTP_STATUS=;
%let HAPLOREG_LD_QUERY_OK=0;
%let HAPLOREG_LD_QUERY_SOURCE=NONE;
%if %sysfunc(countw(%superq(outdsd),.))>1 %then %do;
  %let _out_lib=%scan(%superq(outdsd),1,.);
  %let _out_mem=%scan(%superq(outdsd),-1,.);
%end;
%else %do;
  %let _out_lib=WORK;
  %let _out_mem=&outdsd;
%end;

/*
Optional local cache. The cache is a long table with QUERY_SNP, PROXY_SNP,
LD_POPULATION, PROXY_R2, and CACHE_MIN_R2. It can be produced by
extract_haploreg_ld_for_candidates.pl from the official downloadable archives.
Archive data are available only from r2=0.2 upward; a cache is never used for a
lower requested cutoff. A cache miss falls back to the web query by default.
*/
%let _cache_action=%upcase(%sysfunc(strip(%superq(HAPLOREG_LD_CACHE_MISS_ACTION))));
%if %sysevalf(%superq(_cache_action)=,boolean) %then %let _cache_action=WEB;
%let _cache_min=%sysfunc(inputn(%superq(HAPLOREG_LD_CACHE_MIN_R2),best32.));
%if %sysevalf(%superq(_cache_min)=,boolean) %then %let _cache_min=0.2;

%if not %sysevalf(%superq(HAPLOREG_LD_CACHE_DSD)=,boolean)
    and %sysfunc(exist(%superq(HAPLOREG_LD_CACHE_DSD))) %then %do;
  %if %sysevalf(&ldThresh >= &_cache_min) %then %do;
    data &outdsd;
      length query_snp_rsid rsID $128 chr $16 pos_hg38 r2 8;
      set &HAPLOREG_LD_CACHE_DSD;
      where query_snp=upcase("&snp")
        and ld_population=upcase("&LDpop")
        and proxy_r2 >= &ldThresh;
      query_snp_rsid="&snp";
      rsID=strip(proxy_snp);
      r2=proxy_r2;
      chr='';
      pos_hg38=.;
      keep query_snp_rsid rsID chr pos_hg38 r2;
    run;
    proc sql noprint;
      select count(*) into :_cache_self_n trimmed
      from &outdsd
      where upcase(strip(rsID))=upcase("&snp");
    quit;
    %if %sysevalf(%superq(_cache_self_n)=,boolean) %then %let _cache_self_n=0;
    %if &_cache_self_n>0 %then %do;
      %let HAPLOREG_LD_HTTP_STATUS=LOCAL_CACHE;
      %let HAPLOREG_LD_QUERY_OK=1;
      %let HAPLOREG_LD_QUERY_SOURCE=LOCAL_CACHE;
      %put NOTE: HaploReg LD cache hit snp=&snp population=&LDpop threshold=&ldThresh cache=&HAPLOREG_LD_CACHE_DSD.;
      %return;
    %end;
    %put WARNING: HaploReg LD cache miss for &snp population=&LDpop.;
    %if &_cache_action ne WEB %then %do;
      %let HAPLOREG_LD_HTTP_STATUS=CACHE_MISS;
      %let HAPLOREG_LD_QUERY_SOURCE=LOCAL_CACHE_MISS;
      %return;
    %end;
  %end;
  %else %put WARNING: Ignoring HaploReg LD cache because requested r2=&ldThresh is below cache minimum r2=&_cache_min.;
%end;

%let query_str=%nrstr(&query=)&snp%nrstr(&ldThresh=)&ldThresh%nrstr(&ldPop=)&LDpop%nrstr(&output=text&submit=submit);
%let HAPLOREG_LD_QUERY_SOURCE=WEB;

filename hap temp lrecl=1000000000;
proc http
  url="https://pubs.broadinstitute.org/mammals/haploreg/haploreg.php/post"
  method="POST"
  in="&query_str"
  out=hap
  ct='application/x-www-form-urlencoded';
run;

%let _http_status=&SYS_PROCHTTP_STATUS_CODE;
%let HAPLOREG_LD_HTTP_STATUS=&_http_status;

%if %sysevalf(%superq(_http_status)=,boolean)
    or %sysevalf(&_http_status < 200)
    or %sysevalf(&_http_status >= 400) %then %do;
  %put WARNING: HaploReg LD query failed for &snp population=&LDpop HTTP=&_http_status.;
  data &outdsd;
    length query_snp_rsid rsID $128 chr $16 pos_hg38 r2 8;
    stop;
  run;
  filename hap clear;
  %return;
%end;

%let _response_bytes=0;
data _null_;
  length _size_text $64;
  _fid=fopen('hap');
  if _fid>0 then do;
    _size_text=finfo(_fid,'File Size (bytes)');
    call symputx('_response_bytes',input(compress(_size_text,,'kd'),best32.),'L');
    _rc=fclose(_fid);
  end;
run;
%if %sysevalf(%superq(_response_bytes)=,boolean) %then %let _response_bytes=0;
%if %sysevalf(&_response_bytes <= 0) %then %do;
  %put WARNING: HaploReg returned an empty response for &snp population=&LDpop.;
  data &outdsd;
    length query_snp_rsid rsID $128 chr $16 pos_hg38 r2 8;
    stop;
  run;
  filename hap clear;
  %return;
%end;

/*
PROC IMPORT writes a fatal-looking ERROR when HaploReg returns only a blank
line or a header without data.  Count nonblank records first so an expected
"no LD response" is represented by an empty data set and audit status instead
of poisoning an otherwise successful SAS ODA job log.
*/
%let _response_records=0;
data _null_;
  infile hap lrecl=32767 truncover end=_eof;
  input;
  if not missing(_infile_) then _n_nonblank+1;
  if _eof then call symputx('_response_records',_n_nonblank,'L');
run;
%if %sysevalf(%superq(_response_records)=,boolean) %then %let _response_records=0;
%if %sysevalf(&_response_records < 2) %then %do;
  %put WARNING: HaploReg returned no tabular LD rows for &snp population=&LDpop.;
  data &outdsd;
    length query_snp_rsid rsID $128 chr $16 pos_hg38 r2 8;
    stop;
  run;
  filename hap clear;
  %return;
%end;

proc import datafile=hap dbms=tab out=&outdsd replace;
  getnames=yes;
  guessingrows=max;
run;
filename hap clear;

%if not %sysfunc(exist(&outdsd)) %then %do;
  %put WARNING: HaploReg returned no importable table for &snp population=&LDpop.;
  data &outdsd;
    length query_snp_rsid rsID $128 chr $16 pos_hg38 r2 8;
    stop;
  run;
  %return;
%end;

proc sql noprint;
  select sum(upcase(name)='RSID'),
         sum(upcase(name)='R2')
    into :_has_rsid trimmed,
         :_has_r2 trimmed
  from dictionary.columns
  where libname=upcase("&_out_lib")
    and memname=upcase("&_out_mem")
  ;
quit;

%if %sysevalf(%superq(_has_rsid)=,boolean) %then %let _has_rsid=0;
%if %sysevalf(%superq(_has_r2)=,boolean) %then %let _has_r2=0;

%if &_has_rsid=0 or &_has_r2=0 %then %do;
  %put WARNING: HaploReg response for &snp population=&LDpop lacks rsID or r2 columns.;
  data &outdsd;
    length query_snp_rsid rsID $128 chr $16 pos_hg38 r2 8;
    stop;
  run;
  %return;
%end;

data &outdsd;
  length query_snp_rsid $128;
  set &outdsd(drop=query_snp_rsid);
  query_snp_rsid="&snp";
run;

proc sql noprint;
  select count(*) into :_n_rows trimmed
  from &outdsd
  where upcase(strip(rsID))=upcase("&snp")
  ;
quit;
%if %sysevalf(%superq(_n_rows)=,boolean) %then %let _n_rows=0;
%if %sysevalf(&_n_rows > 0) %then %let HAPLOREG_LD_QUERY_OK=1;

%put NOTE: HaploReg LD query snp=&snp population=&LDpop threshold=&ldThresh HTTP=&_http_status query_ok=&HAPLOREG_LD_QUERY_OK.;
%mend;
