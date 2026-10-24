package BookTrackSVG;
use strict;use warnings;use XML::LibXML;use Exporter 'import';our@EXPORT_OK=qw(inspect_track inspect_cnv);
sub nums{($_[0]// '')=~/[-+]?(?:\d*\.)?\d+(?:e[-+]?\d+)?/ig}
sub mul{my($l,$r)=@_;[$l->[0]*$r->[0]+$l->[2]*$r->[1],$l->[1]*$r->[0]+$l->[3]*$r->[1],$l->[0]*$r->[2]+$l->[2]*$r->[3],$l->[1]*$r->[2]+$l->[3]*$r->[3],$l->[0]*$r->[4]+$l->[2]*$r->[5]+$l->[4],$l->[1]*$r->[4]+$l->[3]*$r->[5]+$l->[5]]}
sub transform{my$s=shift;my$m=[1,0,0,1,0,0];while($s=~/(matrix|translate|scale|rotate)\s*\(([^)]+)\)/g){my($type,$v)=($1,$2);my@v=nums($v);my$t;if($type eq'matrix'){$t=\@v}elsif($type eq'translate'){$t=[1,0,0,1,$v[0],$v[1]//0]}elsif($type eq'scale'){$t=[$v[0],0,0,$v[1]//$v[0],0,0]}else{my$a=$v[0]*atan2(1,1)/45;my$c=cos($a);my$sn=sin($a);$t=[$c,$sn,-$sn,$c,0,0];if(@v==3){$t=mul([1,0,0,1,$v[1],$v[2]],mul($t,[1,0,0,1,-$v[1],-$v[2]]))}}$m=mul($m,$t)}$m}
sub world{my($n,$x,$y)=@_;my$m=[1,0,0,1,0,0];while($n&&$n->nodeType==1){$m=mul(transform($n->getAttribute('transform')//''),$m);$n=$n->parentNode}return[$m->[0]*$x+$m->[2]*$y+$m->[4],$m->[1]*$x+$m->[3]*$y+$m->[5]]}
sub min{my$v=shift;for(@_){$v=$_ if$_<$v}$v}sub max{my$v=shift;for(@_){$v=$_ if$_>$v}$v}
sub inspect_track{
 my($file,$labels)=@_;my$x=XML::LibXML->load_xml(location=>$file);
 my@wall=grep{($_->getAttribute('style')//'')=~/fill: #FFFFFF/&&($_->getAttribute('width')//0)>250&&($_->getAttribute('height')//0)>150&&($_->getAttribute('height')//0)<600}$x->findnodes('//*[local-name()="rect"]');
 die"No bounded plot wall $file"unless@wall;my$w=$wall[-1];my$lt=world($w,$w->getAttribute('x'),$w->getAttribute('y'));my$rb=world($w,$w->getAttribute('x')+$w->getAttribute('width'),$w->getAttribute('y')+$w->getAttribute('height'));
 my(@refs,@guides,@bodies,@circles);
 for my$p($x->findnodes('//*[local-name()="path"]')){my$s=$p->getAttribute('style')//'';my@n=nums($p->getAttribute('d'));next unless@n>=4;my$a=world($p,@n[0,1]);my$b=world($p,@n[2,3]);if($s=~/#CCCCCC/&&@n==4&&abs($a->[1]-$b->[1])<.001&&abs($a->[0]-$b->[0])>.8*($rb->[0]-$lt->[0])){push@refs,$a->[1]}
 if($s=~/stroke-dasharray/&&@n==4&&abs($a->[0]-$b->[0])<.001&&abs($a->[1]-$b->[1])>50){push@guides,{x=>$a->[0],top=>min($a->[1],$b->[1]),bottom=>max($a->[1],$b->[1])}}
 my($width)=$s=~/stroke-width:\s*([0-9.]+)/;if(defined($width)&&$width>=7&&$s=~/stroke-opacity/){my$u=world($p,0,0);my$v=world($p,0,$width);push@bodies,{x1=>$a->[0],x2=>$b->[0],y=>$a->[1],thickness_px=>abs($v->[1]-$u->[1])}}
 }
 for my$c($x->findnodes('//*[local-name()="circle"]')){next unless($c->getAttribute('style')//'')=~/#(?:0072B2|D55E00|CC79A7|009E73)/i;my$p=world($c,$c->getAttribute('cx'),$c->getAttribute('cy'));next unless$p->[0]>=$lt->[0]&&$p->[0]<=$rb->[0]&&$p->[1]>=$lt->[1]&&$p->[1]<=$rb->[1];my$edge=world($c,$c->getAttribute('cx')+$c->getAttribute('r'),$c->getAttribute('cy'));push@circles,{x=>$p->[0],y=>$p->[1],radius_px=>abs($edge->[0]-$p->[0])}}
 my$zero=max(@refs);my$top=min(@refs);@bodies=grep{$_->{y}>$zero}@bodies;my@l;
 for my$label(@{$labels//[]}){my@n=grep{$_->textContent eq$label}$x->findnodes('//*[local-name()="text"]');die"Expected unique $label in $file"unless@n==1;my$t=$n[0];my$len=$t->getAttribute('textLength');die"Missing emitted length $label"unless length($len//'');my$fs=$t->getAttribute('style')//'';my($font)=$fs=~/font-size:\s*([0-9.]+)/;$font//=11;my($tx,$ty)=map{$t->getAttribute($_)||0}qw(x y);my@pts=map{world($t,@$_)}([$tx,$ty-$font*.92],[$tx+$len,$ty-$font*.92],[$tx,$ty+$font*.25],[$tx+$len,$ty+$font*.25]);my$bounds=[min(map{$_->[0]}@pts),min(map{$_->[1]}@pts),max(map{$_->[0]}@pts),max(map{$_->[1]}@pts)];push@l,{label=>$label,bounds_px=>$bounds,above_wall_px=>$bounds->[1]-$lt->[1],above_separator_px=>$top-$bounds->[3],center_error_px=>($bounds->[1]+$bounds->[3]-$lt->[1]-$top)/2};
 }
 my@names;for my$t($x->findnodes('//*[local-name()="text"]')){my$s=$t->getAttribute('style')//'';next unless$s=~/font-style: italic/;my($font)=$s=~/font-size:\s*([0-9.]+)/;$font//=11;my($tx,$ty)=map{$t->getAttribute($_)||0}qw(x y);my$len=$t->getAttribute('textLength')||0;my@pts=map{world($t,@$_)}([$tx,$ty-$font*.92],[$tx+$len,$ty-$font*.92],[$tx,$ty+$font*.25],[$tx+$len,$ty+$font*.25]);my$bounds=[min(map{$_->[0]}@pts),min(map{$_->[1]}@pts),max(map{$_->[0]}@pts),max(map{$_->[1]}@pts)];my@near=grep{$_->{y}>=$bounds->[3]&&$bounds->[0]<=$_->{x2}+1&&$bounds->[0]>=$_->{x1}-1}@bodies;my@clear=map{$_->{y}-$_->{thickness_px}/2-$bounds->[3]}@near;push@names,{label=>$t->textContent,bounds_px=>$bounds,body_clearance_px=>@clear?min(@clear):undef,separator_clearance_px=>$bounds->[1]-$zero}}
 my@alignment;for my$g(@guides){my@d=map{abs($_->{x}-$g->{x})}@circles;push@alignment,{guide_x=>$g->{x},nearest_dot_x_residual_px=>@d?min(@d):undef}}
 my@view=nums($x->documentElement->getAttribute('viewBox'));
 return{file=>$file,svg_viewbox=>\@view,wall_px=>[@$lt,@$rb],headroom_separator_y_px=>$top,zero_separator_y_px=>$zero,annotation_share=>($rb->[1]-$zero)/($rb->[1]-$lt->[1]),guides=>\@alignment,marker_count=>scalar@circles,markers=>\@circles,body_segments=>\@bodies,gene_labels=>\@names,labels=>\@l};
}
sub inspect_cnv {
 my($file)=@_;my$x=XML::LibXML->load_xml(location=>$file);my(@circles,@intervals);
 for my$c($x->findnodes('//*[local-name()="circle"]')){my$p=world($c,$c->getAttribute('cx'),$c->getAttribute('cy'));push@circles,{x=>$p->[0],y=>$p->[1],style=>$c->getAttribute('style'),radius=>$c->getAttribute('r')}}
 for my$p($x->findnodes('//*[local-name()="path"]')){my$s=$p->getAttribute('style')//'';next unless$s=~/stroke: #FF8C00;/;my@n=nums($p->getAttribute('d'));die'Unexpected CNV path'unless@n==4;my$a=world($p,@n[0,1]);my$b=world($p,@n[2,3]);push@intervals,{x1=>$a->[0],x2=>$b->[0],y=>$a->[1],style=>$s}}
 return{file=>$file,circles=>\@circles,intervals=>\@intervals};
}
1;
