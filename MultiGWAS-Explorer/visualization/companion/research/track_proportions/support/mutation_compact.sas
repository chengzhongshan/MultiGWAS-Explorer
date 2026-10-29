/* Named book-owned display extension after the original mutation wrapper.
   Input must be the checked output of book_track_display, not raw mutations. */
%macro BookMutationTrackDisplay(data=work.track_display,
 out=work.mutation_compact_display,height=350,dots=11,label_lift=0.30);
data &out;
 set &data;
 /* Keep the gene name clear of the enlarged body stroke. */
 if not missing(book_gene_y) and book_gene_y<0 then
  book_gene_y=book_gene_y+&label_lift;
run;
proc template;source BedGraph / file="&BOOK_OUT/mutation_compact_resolved.sas";run;
data _null_;
 infile "&BOOK_OUT/mutation_compact_resolved.sas" lrecl=32767 truncover;
 file "&BOOK_OUT/mutation_compact_template.sas" lrecl=32767;
 length line $32767;input;line=_infile_;
 line=prxchange("s/designheight=[0-9.]+/designheight=&height/",1,line);
 if index(lowcase(line),'drawtext textattrs=') and index(lowcase(line),'rotate=90') then
  line=tranwrd(line,'rotate=90;','rotate=90 width=100 widthunit=pixel justify=center;');
 if index(lowcase(line),'scatterplot x=book_marker_pos ') then
  line=prxchange("s/size=[0-9.]+(?:px|pt)?/size=&dots.px/",1,line);
 put line;
run;
proc template;%include "&BOOK_OUT/mutation_compact_template.sas";run;
%mend;
