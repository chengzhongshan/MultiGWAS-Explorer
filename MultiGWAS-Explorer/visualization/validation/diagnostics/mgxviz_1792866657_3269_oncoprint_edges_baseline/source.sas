/* Invented 2x2 assay: duplicate event, multiple classes, 100% and 0/0. */
data work.samples;length sample $3 cohort $8;
 input sample $ sample_index cohort $;datalines;
C1 1 ArmA
C2 2 ArmB
;
run;
data work.genes;length gene $8;input gene $ gene_index;datalines;
GENE1 1
GENE2 2
;
run;
data work.calls;length sample $3;input sample $ gene_index callable;datalines;
C1 1 1
C2 1 1
C1 2 0
;
run;
data work.variants;length sample $3 gene $8 variant_id $10 mutation_class $12;
 input sample $ gene $ variant_id $ mutation_class $;datalines;
C1 GENE1 e1 Missense
C1 GENE1 e1 Missense
C1 GENE1 e2 Truncating
C2 GENE1 e3 Splice
;
run;
%MGXOncoPrint(samples=work.samples,genes=work.genes,
 callability=work.calls,variants=work.variants,
 prefix=edge,outdir=&BOOK_OUT,image=oncoprint_edges,width=6,height=4,
 template_file="oncoprint_template.sas");
data work.edge_status;status=&MGX_ONCOPRINT_RC;run;
