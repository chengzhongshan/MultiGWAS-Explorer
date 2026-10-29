%macro make_fake_axis_values4grps(
/*This macro has issue when axis_var containing positve and negative values among different grps*/
dsdin,
axis_var,
axis_grp,
new_fake_axis_var,
dsdout,
yaxis_macro_labels=ylabelsmacro_var,
step4yaxis_macro_labels=1,/*Only keep ticks with mod(t,step)=0, which will prevent the y-axis ticks from too compacted with each other!*/
fc2scale_pos_vals=2, /*Use this fc to enlarge the proportion of positive values in the plots
It seems that fc=2 is the best for the final ticks of different tracks;
*/
mod_num2keep= /*Default is empty for not filtering these elements by mod; when values, 
such as 2 or 3 are provided, only keep numbers that fulfil the mod(element,num)=0*/
);

*Fix a bug when only one record in a group;
*make an extra copy for the records to make it has >1 records;
*This is because the script will not be able to update corrected;
*y axis labels if there is only one record for a group;
proc sql;
create table _single_ as 
select *
from &dsdin
group by &axis_grp
having count(&axis_grp)=1;
data &dsdin;
set &dsdin _single_;
data &dsdin;
set &dsdin;
*Scale position value and make it larger by fold change;
if &axis_var>0 then &axis_var=&fc2scale_pos_vals*&axis_var;
run;

*add 2 to the max value to make grps separated better;
proc sql;
create table &dsdout as
select a.*,&axis_var as &new_fake_axis_var,
       floor(min(&axis_var)) as min_val,
       ceil(max(&axis_var))+2 as max_val
from &dsdin as a
group by &axis_grp
order by &axis_grp,&axis_var;

data &dsdout;
retain max grpnum mtag 0;
set &dsdout;
grp_end_tag=0;
if first.&axis_grp then do;
*important to set max^=0;
*No other value can be used, as 0 is very special;
 if max^=0 then do;
/*	 if min_val<0 then do;*/
/*		 *The 1st one of each grp would be the min_val;*/
/*   &new_fake_axis_var=max-min_val;*/
/*		end;*/
/*		else do;*/
/*   &new_fake_axis_var=&axis_var+max;*/
/*		end;*/
	  *As max value is integer, the 1st y value if not integer, need to add back its decimal value;
	  &new_fake_axis_var=max+(&axis_var-min_val);
   grpnum=grpnum+1;
 end;
 else do;
  &new_fake_axis_var=&axis_var;
		mtag=1;
  grpnum=1;
 end;
end;

else if last.&axis_grp then do;
 grp_end_tag=1;
	*important to set max^=0;
	*also need to check mtag,which forces to run the following when passed the;
	*the 1st element;
 if max^=0 then do;
   &new_fake_axis_var=(&axis_var-min_val)+max;
/*		*for negative min_val, it is necessary to adjust max;*/
/*		if min_val<0 then do;*/
/*		  *when max is still equal to 0, keep it as 0;*/
/*		  if max^=0 then do;*/
/*     max=max+max_val-min_val;*/
/*					*Important to make it as the close largest integer;*/
/*					max=ceil(max);*/
/*				end;*/
/*		end;*/
/*		else do;*/
/*			 *max=max+max_val;*/
/*		  *The above can not be used;*/
/*		  *Not 100% sure why the following works!;*/
/*		   max=&axis_var-min_val+2+max;*/
/*					max=ceil(max);*/
/*		end;*/
			 *Use correct syntax;
			 max=ceil(&new_fake_axis_var)+2;
 end;
 else do;
  &new_fake_axis_var=&axis_var;
		*for negative min_val, it is necessary to adjust max;
		*Only if max^=0, the value of max will be updated;
		if min_val<0 and max^=0 then do;
   max=max+max_val-min_val;
			max=ceil(max);
		end;
		else do;
			max=max+max_val;
			max=ceil(max);
		end;
 end;
	mtag=0;
end;

else do;
  if min_val<0 then do;
     &new_fake_axis_var=&axis_var-min_val+max;
		end;
		else do;
		  *this is wrong;
/*				&new_fake_axis_var=&axis_var+max;*/
	  	&new_fake_axis_var=&axis_var-min_val+max;
		end;																																																	                                                                                                                                                           
		if max=0 and mtag=1 then &new_fake_axis_var=&axis_var;
end;
output;
*drop min max min_val max_val;
*Make sure to sort by two vars;
by &axis_grp &axis_var;
run;

*Generate real y axis labels;
%global &yaxis_macro_labels fake_max_y fake_min_y fake_refline_values;
proc sql noprint;
select max_val,min_val,max-1 
into: max_y4grps separated by " ",
    : min_y4grps separated by " ",
				: fake_refline_values separated by " "
from &dsdout
where grp_end_tag=1
order by grpnum,&axis_var;

%put max_y4grps: &max_y4grps;
%put min_y4grps: &min_y4grps;

*Also get the min_val when grp_end_tag=1;
*Get the macro var fake_max_y value;
data _null_;
set &dsdout end=eof;
if eof then call symputx('fake_max_y',ceil(&new_fake_axis_var));
run;

%let fake_min_y=%scan(&min_y4grps,1,%str( ));
%let fake_min_y=%sysfunc(floor(&fake_min_y));

%let nums=;
%do xi=1 %to %sysfunc(countw(&max_y4grps));
  *use %str( ) to prevent from lossing of negative nums;
		%let _min_y=%scan(&min_y4grps,&xi,%str( ));
		%put _min_y is &_min_y;
		%if %eval(&xi=1) %then %do;
   /* %nums_in_range(st=&_min_y,end=%eval(%scan(&max_y4grps,&xi,%str( ))-1),by=1,outmacrovar=nums&xi,quote=1);*/
    %nums_in_range_adj_scale(st=&_min_y,end=%eval(%scan(&max_y4grps,&xi,%str( ))-1),by=1,outmacrovar=nums&xi,
                   filter4scaledvals=%str(>0),scale=&fc2scale_pos_vals,quote=1,mod_num2keep=&mod_num2keep);
				%let nums=&&nums&xi;
				*Need to replace the last value as empty;
   %let nums=%sysfunc(prxchange(s/\S+$/" "/,-1,&nums));
		%end;
		%else %do;
/*		  %nums_in_range(st=&_min_y,end=%eval(%scan(&max_y4grps,&xi,%str( ))-1),by=1,outmacrovar=nums&xi,quote=1); */
		  %nums_in_range_adj_scale(st=&_min_y,end=%eval(%scan(&max_y4grps,&xi,%str( ))-1),by=1,outmacrovar=nums&xi,
                                 filter4scaledvals=%str(>0),scale=&fc2scale_pos_vals,quote=1,mod_num2keep=&mod_num2keep);
				*Need to replace the last value as empty;
				%let new_nums=%sysfunc(prxchange(s/\S+$/" "/,-1,&&nums&xi));
				%put modified new_nums are: &new_nums;
				%let nums=&nums &new_nums;
		%end;
	%end;
	%let &yaxis_macro_labels=&nums;
	%put generated the global macro var &&yaxis_macro_labels for labeling y axis, which are: &&&yaxis_macro_labels;
		%put generated the global macro var fake_min_y, the value of which is &fake_min_y;
	 %put generated the global macro var fake_max_y, the value of which is &fake_max_y;
		%put generated the global macro var fake_refline_values, the value of which can be used to make reflines to separate grps:;
		%put &fake_refline_values;
%mend;
