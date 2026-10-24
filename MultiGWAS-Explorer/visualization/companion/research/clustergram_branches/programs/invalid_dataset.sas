/* Deliberate failure: no plot should be produced. The original macro checks
   dataset existence before attempting to traverse a tree. */
%subbranches(data=work.no_such_tree,
  y_name_var=x_name_,y_parent_var=x_parent_,y_height_var=x_height_,
  out=work.should_not_exist,parent=CL10,render=NO);
