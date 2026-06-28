/* libname FM '/home/cheng.zhong.shan/my_shared_file_links/cheng.zhong.shan/F_vs_M_Covid19_Hosp'; */
/* *hg19 version; */
/* %let gtf_gz_url=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz; */
/* %get_genecode_gtf_data(gtf_gz_url=&gtf_gz_url,outdsd=gtf_hg19); */
/* *Or use the hg38 version; */
/* %let gtf_gz_url=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz; */
/* %get_genecode_gtf_data(gtf_gz_url=&gtf_gz_url,outdsd=gtf_hg38); */

/* proc datasets nolist; */
/* copy in=work out=FM memtype=data move; */
/* *select gtf_hg19; */
/* select gtf_hg38; */
/* run; */
/*%macroparas(macrorgx=zip);*/
x cd "G:\NGS_lib\Linux_codes_SAM\Conda_and_Docker_Related_Scripts\perlMCP4Gemini_Paper";

%include "Manhattan4DiffGWASs_png.sas";
/*%debug_macro;*/
%ImportFileHeadersFromZIP( 
zip=E:\LongCOVID_HGI_GWAS\PGC_Large_GWASs\PGC_SCZ_Sex_Stratified_GWASs\scz_mh_p005_chr23.tsv.gz,/*Only provide file with .gz, .zip, or common text file without comporession 
Note: it is necessary to have fullpath for the input file!*/ 
filename_rgx=., 
obs=max, 
sasdsdout=scz_mh, 
deleteZIP=0, 
infile_command=%str( 
dlm='09'x dsd firstobs=2 truncover lrecl=32767;
input CHR
    BP
    SNP :$20.
    ALL_STD_P
    ASN_STD_P
    EUR_STD_P
    ALL_DIFF_P
    ASN_DIFF_P
    EUR_DIFF_P;), 
/*Better to use nrbquote to replace str and use unquote within the macro 
to get back the input infile_command;*/ 
extra_infile_macrovar_prefix=infile_cmd,/*To prevent the crash of sas when the length of the macro var infile_command is too long, 
it is better to assign different parts of infile commands into multiple global macro vars with similar prefix, such as infile_cmd; 
it is better to use bquote or nrbquote to excape each extra infile command!*/ 
num_infile_macro_vars=0,/*Provide positve number to work with the global macro var of extra_infile_macrovar_prefix*/ 
use_zcat=0, 
var4endlinenum=adj_endlinenum, /*make global var for the endline number but it is 
necessary to use syminputx in the infile_command to record the endline number; 
call symputx("&var4endlinenum",trim(left(put(_n_,8.)))); 
It is possible to assign other numeric value generated in the infile_command to 
this macro var for other purpose, because this global macro var will be accessible 
by other outsite macros! 
call symputx('adj_endlinenum',trim(left(put(rowtag,8.))));*/ 
global_var_prefix4vars2drop=drop_var,/*To handle the issue of trunction of macro var infile_command if there are too many variables to be dropped in the infile procedure; 
it is feasible to create global macro variables with the same prefix, such as drop_var, to exclude them*/ 
num_vars2drop=0 /*Provide postive number to work with the macro var global_var_prefix4vars2drop to resolve these variables to be excluded*/ 
); 

/*
filename mhgz zip "scz_mh_p005_chr23.tsv.gz" gzip;
data scz_mh;
  infile mhgz dlm='09'x dsd firstobs=2 truncover lrecl=32767;
  length SNP $40;
  input
    CHR
    BP
    SNP :$40.
    ALL_STD_P
    ASN_STD_P
    EUR_STD_P
    ALL_DIFF_P
    ASN_DIFF_P
    EUR_DIFF_P
  ;
run;

proc sort data=scz_mh;
  by CHR BP;
run;
*/

%Manhattan4DiffGWASs(
  dsdin=scz_mh,
  pos_var=BP,
  chr_var=CHR,
  P_var=ALL_STD_P,
  Other_P_vars=ASN_STD_P EUR_STD_P ALL_DIFF_P ASN_DIFF_P EUR_DIFF_P,
  logP=1,
  gwas_thrsd=7.30103,
  dotsize=1,
  _logP_topval=10,
  y_axix_step=2,
  fig_width=1200,
  fig_height=700,
  fontsize=3,
  gwas_label_names=%str(All standardized P|Asian standardized P|European standardized P|All female-vs-male diff P|Asian female-vs-male diff P|European female-vs-male diff P),
  gwas_label_x_pct=50,
  gwas_label_y_frac=0.90,
  gwas_label_size=3.6,
  gwas_label_halo_size=4.9,
  gwas_label_angle=0,
  flip1stGWAS_signal=0,
  rm_signals_with_logP_lt=0.5,
  outputfigname=PGC_SCZ_SAS_manhattan,
  Use_scaled_pos=1,
  sep_chr_grp=0,
  gwas_sortedby_numchrpos=1
);


libname FM "E:\LongCOVID_HGI_GWAS";
data out;
set scz_mh;
where (ASN_STD_P>0 and ASN_STD_P<1e-6) or (EUR_STD_P>0 and EUr_STD_P<1e-6);
run;

/*%macroparas(macrorgx=hap);*/
%get_top_signal_within_dist( 
dsdin=out, 
grp_var=chr, 
signal_var=ASN_STD_P, 
select_smallest_signal=1, 
pos_var=BP, 
pos_dist_thrshd=1e8, 
dsdout=top_hit4diffp, 
signal_thrshd=1e-6 /*filter the input dsdin by association P, i.e, &signal_val <= &signal_thrshd*/ 
); 
*Only focus on SNPs with both nominal significant association with both traits;
/*
data top_hit4diffp;
set top_hit4diffp;
where gwas1_p<=0.05 and gwas2_p<=0.05 and pval<=5e-8;
run;
*/
proc sql noprint;
select SNP into: top_snps separated by ' '
from top_hit4diffp 
order by chr,BP;

%QueryHaploreg(/*Query Haploreg4 for each input SNP to get genes close to it!*/ 
rsids=&top_snps, 
dsdout=snps2genes, 
print_html=1 /*Print out the annotations of query SNP(s)*/ 
);
*Add genes to these top SNPs;
proc sql;
create table top_hit4diffp as
select a.*,b.gene 
from top_hit4diffp as a
left join
snps2genes as b
on a.SNP=b.rsid;
 
create table top_sigs as 
select a.*,catx(':',b.SNP,b.gene) as snp_gene
from scz_mh as a,
         top_hit4diffp as b
where a.chr=b.chr and a.BP between (b.BP-1e7) and (b.BP+1e7);

select catx(":",SNP,gene) into: snp_gene_label separated by ' '
from top_hit4diffp
order by chr,BP;

/*
%ds2csv(data=top_hit4diffp,
csvfile=E:\LongCOVID_HGI_GWAS\PGC_Large_GWASs\PGC_GWAS_Analyzer_Paper\DiffGWAS4PTSD_vs_SCZ_tophits2genes.csv,
runmode=b);

data G.top_sigs4PGC;
set top_sigs;
run;
*/

%Manhattan4DiffGWASs( 
dsdin=top_sigs, 
pos_var=BP, 
chr_var=snp_gene, 
P_var=ALL_STD_P,
Other_P_vars=ASN_STD_P EUR_STD_P ALL_DIFF_P ASN_DIFF_P EUR_DIFF_P,
gwas_thrsd=7.3,/*Use it to draw significance reference line in each GWAS track*/ 
thrsd_line_color=gray,
dotsize=1,/*The dot size for scatter plots*/
_logP_topval=8, /*Top -log10P value to truncate GWAS signals and also restrict the max yaxis value of each GWAS track;
Make sure to input EVEN number for the macro, as the macro separate ticks by step 2!*/
y_axix_step=5,/*Customize the step for all y-axis tikets*/
rm_signals_with_logP_lt=0, 
flip1stGWAS_signal=0, 
fig_width=1200, 
fig_height=1000, 
angle4xaxis_label=90, 
xgrp_y_pos=-5, 
yoffset_setting=%str(offset=(15,0.5)),
draw_local_Manhattan=0, 
sep_chr_grp=1, /*Default is not to add lines to separate x-axis chromosomal groups;
target_SNPs=&top_snps,/*Default is empty; please provide rsid that can be matched with the snp macro variable*/
snp_var=rsid,/*It is necessary to have snp_var supplied when drawing local Manhattan plot*/
snp_gene_splitter=:,/*In case the gene name for the snp is also supplied to the snp_var, the macro
will split the snp_var into two string, with the first is snp and the 2nd is genename for it, which will
be plotted at the bottom of the figure as x-axis labels for different snp mahattan plot*/
target_SNPs=&snp_gene_label,/*Default is empty; please provide rsid that can be matched with the snp macro variable*/
Keep_order_of_target_SNPs=1, /*Draw local Manhattan plot according to the order of target SNPs
Note: need to set this macro with value 1 if drawing local Manhattan plots for target SNPs or top hits, 
which means if either target_SNPs or top_hit_thresd is not empty, please assign value 1 to this macro var!
When draw genome-wide Manhattan plots, it is required to assign value 1 to this macro var.*/
top_hit_thresd=1e-6,/*provide a p value threshold to only draw local Manhattan plot for the smallest 
top hit around a specific genomic window,such as p < 1e-6 within a window of 1e7 bp*/
dist4get_smallest_top_hit=1e7, /*Select the smallest top SNP around a genomic window of the supplied distance in bp*/
only_get_top_hit4n_th_gwas=3 /*The parameter enables the macro to focus on top hits from specific gwas represented by 
its order starting from 1 to n for the 1st gwas and other gwass inferred by their supplied P variables; the default value 0 means
to query all gwas top hits; if only want to query top hits from the 1 gwas, please supply value 1, and this applicable to 
other gwass if the correct numeric order for the gwas is supplied here!*/ 
);

***Important:
*The following codes need to be modifed to use the corresponding variables in the input dsd for the macro Manhattan4DiffGWASs, which will be used to draw the local Manhattan plot for top hits;
*Note: SNP_Local_Manhattan_With_GTF is accessible within SAS ODA automatically;

*Run analysis here based on previously saved results;
%let all_snps=&top_snps;
data top_sigs;
set top_sigs;
*Get Z-scores for both gwass for comparable coloring of scatter plots in local Manhattan plot;
gwas1_z=gwas1_beta/gwas1_se;
gwas2_z=gwas2_beta/gwas2_se;
run;
%SNP_Local_Manhattan_With_GTF(/*As this macro use other sub-macros, it is not uncommon that some global macro
vars would be in the same name, such as macro vars chr and i, thus, to avoid of crash, chr_var is used instead of macro
var chr in this macro*/
gwas_dsd=top_sigs,
chr_var=chr,
AssocPVars=%pull_list(input_list=gwas1_p gwas2_p pval,idx4pull=1 2 3),
SNP_IDs=&all_snps,
/*if providing chr:pos or chr:st:end, it will query by pos;
Please also enlarge the dist2snp to extract the whole gene body and its exons,
altought the final plots will be only restricted by the input st and end positions!*/
dist2snp=500000,
/*in bp; left or right size distant to each target SNP for the Manhattan plot*/
SNP_Var=rsid,
Pos_Var=pos,
gtf_dsd=FM.GTF_HG38,/*Reuse the gtf data in SAS ODA by searching it*/,
ZscoreVars=%pull_list(input_list=gwas1_z gwas2_z diff_zscore,idx4pull=1 2 3),/*Can be beta1 beat2 or other numberic vars indicating assoc or other +/- directions*/ 
gwas_labels_in_order=%pull_list(input_list=PGC_Schizophrenia PGC_PTSD Schizophrenia_vs_PTSD,idx4pull=1 2 3),/*If providing _ for labeling each GWAS, 
the _ will be replaced with empty string, which is useful when wanting to remove gwas label 
if only one scatterplot or the label for a gwas containing spaces;
The list will be used to label scatterplots 
by the sub-macro map_grp_assoc2gene4covidsexgwas*/
design_width=950, /*475*600 is the best width and height for publication*/
design_height=800, 
barthickness=10, /*gene track bar thinkness*/
dotsize=5, 
dist2sep_genes=20000000,/*Distance to separate close genes into different rows in the gene track; provide negative value
to have all genes in a single row in the final gene track*/
where_cndtn_for_gwasdsd=%str(), /*where condition to filter input gwas_dsd*/

shift_text_yval=0.2, /*in terms of gene track labels, add positive or negative vale, ranging from 0 to 1, 
                      to liftup or lower text labels on the y axis; the default value is -0.2 to put gene lable under gene tracks;
                      Change it with the macro var pct4neg_y!*/
fig_fmt=png, /*output figure formats: svg, png, jpg, and others*/
pct4neg_y=2, /*the most often used value is 1;
              compacting the bed track y values by increasing the scatterplot scale, 
              which can reduce the bed trace spaces; It seems that two-fold increasement
              leads to better ticks for different tracks!
              Use value >1 will increase the gene tract, while value < 1 will reduce it!
              Note: when there are only 1 or 2 scatterplots, it is better to provide value = 0.5;
              Modify this parameter with the parameter shift_text_yval to adjust gene label!
              Typically, when there are more scatterplots, it is necessary to increase the value of pct4neg_y accordingly;
              If there are only <4 scatterplots, the value would be usually set as 1 or 2;
              */
adjval4header=-2, /*In terms of header of each subscatterplot, provide postive value to move up scatter group header by the input value*/
makedotheatmap=1,/*use colormap to draw dots in scatterplot instead of the discretemap;
Note: if makedotheatmap=1, the scatterplot will not use the discretemap mode based on
the negative and postive values of lattice_subgrp_var to color dots in scatterplot*/

color_resp_var=,/*Use value of the var to draw colormap of dots in scatterplot
if empty, the default var would be the same as that of yval_var;*/

makeheatmapdotintooneline=0,/*This will make all dots have the same yaxis value but have different colors 
based on its real value in the heatmap plot; To keep the original dot y axis value, assign 0 to the macro var
This would be handy when there are multiple subgrps represented by different y-axis values! By modifying
the y-axis values for these subgrps, the macro can plot them separately in each subtrack!
*/
var4label_scatterplot_dots= ,/*Make sure the variable name is not grp, which is a fixed var used by the macro for other purpose;
the variable should contain values of target SNPs and other non-targets are asigned with empty values;
Whenever  makeheatmapdotintooneline=1 or 0, it is possible to use values of the var4label_scatterplot_dots to
label specific scatterplot dots based on the customization of the variable predifined by users for the input data set; 
default is empty; provide a variable that include non-empty strings for specific dots in the 
scatterplots;*/
SNPs2label_scatterplot_dots=&all_snps, /*Add multiple SNP rsids to label dots within or at the top of scatterplot
Note: if this parameter is provided, it will replace the parameter var4label_scatterplot_dots!
*/
yoffset4max_drawmarkersontop=0.3, /*If draw scatterplot marker labels on the top of track, 
 this fixed value will be used instead of yaxis_offset4max!*/
Yoffset4textlabels=5, /*Move up the text labels for target SNPs in specific fold*/
verbose=0
);
