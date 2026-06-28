/*
%Boxplots4GenesInGTExV8ByGrps(
genes=EXOC6B CYP26B1 DYSF,
GTEx_AA_EA_info=H:\F_Queens\360yunpan\SASCodesLibrary\SAS-Useful-Codes\DemoCodes4Macros\GTEX_AA_EA_others.txt,
dsdout=exp,
bygrps=sex AA, 
UseGeneratedDsd=0,
PreviousDsd=tgt,
Lib4PreviousDsd=GTEx,
WhereFilters4Boxplot=%str(),
boxplot_width=1200,
boxplot_height=1000,
draw_exp_heatmap=1, 
columns=1, 
uniscale_type=all, 
yaxis_max_value= 
);
*/

*Go to the SAS work directory and copy all data into the following directory;
*E:\LongCOVID_HGI_GWAS\GTEx_Exp_SAS_Data;
*Run codes here based on previously generated data;
/*%debug_macro;*/
 libname GTEx "E:\LongCOVID_HGI_GWAS\GTEx_Exp_SAS_Data";
 %let genes=EXOC6B CYP26B1 DYSF;
/* %let genes=CRHR1 MAPT;*/
%Boxplots4GenesInGTExV8ByGrps(
genes=&genes,
GTEx_AA_EA_info=H:\F_Queens\360yunpan\SASCodesLibrary\SAS-Useful-Codes\DemoCodes4Macros\GTEX_AA_EA_others.txt,
dsdout=exp,
bygrps=sex AA, /*by AA population and sex or either of them!*/
UseGeneratedDsd=1,
PreviousDsd=tgt, /*This is a fixed dsd used by the proc sgpanel*/
Lib4PreviousDsd=GTEx,
WhereFilters4Boxplot=%str(prxmatch("/Brain/",cluster) and prxmatch("/(EA)/",AA_Sex);),/*such as cluster in ("Lung" "Spleen" "Whole Blood"); Note these items are case sensitive!*/
boxplot_width=1000,
boxplot_height=1600,
draw_exp_heatmap=1, /*Also draw exp heatmap*/
columns=1, /*Draw boxplots in the number of columns
asign value >1 if drawing multiple genes column-wide!*/
uniscale_type=all, /*make the axis the same for ALL, column (column-wide), or row (row-wide)*/
yaxis_max_value=, /*set the max value for the y axis, which will be set for all lattice panels
If left empty, the macro will get the closer largest integer as the max y value;*/
WhereFilters4Heatmap=%str(prxmatch("/(EA)/",AA_Sex)),/*such as cluster in ("Lung" "Spleen" "Whole Blood") or 
%str(prxmatch("/Brain/",cluster) and prxmatch("/(EA|AA)/",AA_Sex)), without adding ending code ";"
Note: this filters only applicable to heatmaps using proc sql where condition!*/
heatmap_height=35, /*heatmap height in cm*/
heatmap_width=20,/*default is empty, letting the macro to determined it based on the total number of genes*/
heatmapcolumnweights=0.05 0.95, /*figure 2 column ratio*/
heatmaprowweights=0.05 0.95, /*figure 2 row ratio*/
ht_xtick_modifiers=%nrbquote(s/F.EA/Female/) %nrbquote(s/M.EA/Male/)
);
libname GTEx clear;



 libname GTEx "E:\LongCOVID_HGI_GWAS\GTEx_Exp_SAS_Data";
/* %let genes=EXOC6B CYP26B1 DYSF;*/
 %let genes=CRHR1 MAPT;
%Boxplots4GenesInGTExV8ByGrps(
genes=&genes,
GTEx_AA_EA_info=H:\F_Queens\360yunpan\SASCodesLibrary\SAS-Useful-Codes\DemoCodes4Macros\GTEX_AA_EA_others.txt,
dsdout=exp,
bygrps=sex AA, /*by AA population and sex or either of them!*/
UseGeneratedDsd=1,
PreviousDsd=tgt, /*This is a fixed dsd used by the proc sgpanel*/
Lib4PreviousDsd=GTEx,
WhereFilters4Boxplot=%str(prxmatch("/Brain/",cluster) and prxmatch("/(EA)/",AA_Sex);),/*such as cluster in ("Lung" "Spleen" "Whole Blood"); Note these items are case sensitive!*/
boxplot_width=800,
boxplot_height=800,
draw_exp_heatmap=1, /*Also draw exp heatmap*/
columns=1, /*Draw boxplots in the number of columns
asign value >1 if drawing multiple genes column-wide!*/
uniscale_type=all, /*make the axis the same for ALL, column (column-wide), or row (row-wide)*/
yaxis_max_value=, /*set the max value for the y axis, which will be set for all lattice panels
If left empty, the macro will get the closer largest integer as the max y value;*/
WhereFilters4Heatmap=%str(prxmatch("/(EA)/",AA_Sex)),/*such as cluster in ("Lung" "Spleen" "Whole Blood") or 
%str(prxmatch("/Brain/",cluster) and prxmatch("/(EA|AA)/",AA_Sex)), without adding ending code ";"
Note: this filters only applicable to heatmaps using proc sql where condition!*/
heatmap_height=35, /*heatmap height in cm*/
heatmap_width=25,/*default is empty, letting the macro to determined it based on the total number of genes*/
heatmapcolumnweights=0.05 0.95, /*figure 2 column ratio*/
heatmaprowweights=0.05 0.95, /*figure 2 row ratio*/
ht_xtick_modifiers=%nrbquote(s/F.EA/Female/) %nrbquote(s/M.EA/Male/)
);
libname GTEx clear;
