
%macro RandBetween(min, max);
%local min max rnd _unif_num_ _rnd_;
%let _unif_num_=%sysfunc(ranuni(0));
%let _rnd_=%sysevalf((1+&max-&min)*&_unif_num_);
%let rnd=%sysfunc(floor(&_rnd_));
&rnd
%mend;