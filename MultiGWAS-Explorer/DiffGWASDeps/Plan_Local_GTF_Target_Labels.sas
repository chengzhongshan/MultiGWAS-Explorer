/* SAS ODA fallback for the Perl local-GTF label planner. This macro runs on
   the imported locus only; it changes label positions, never variant BPs. */
%macro Plan_Local_GTF_Target_Labels(
  wide_dsd=,
  label_snps=,
  center_bp=,
  window_bp=,
  design_width=950,
  design_height=1000,
  font_size=10,
  requested_layout=auto,
  out_angle=gtf_label_text_rotate_angle,
  out_positions=gtf_label_positions,
  out_headroom=gtf_label_headroom_frac,
  out_center_offset=gtf_label_center_offset
);
  %if %length(%superq(wide_dsd))=0 or %length(%superq(label_snps))=0 %then %do;
    %put ERROR: SAS local-GTF label planner needs a wide dataset and SNP labels.;
    %abort 255;
  %end;
  proc sort data=&wide_dsd(keep=SNP BP where=(not missing(SNP) and not missing(BP)))
    out=_gtf_label_candidates nodupkey;
    by SNP BP;
  run;
  data _gtf_label_candidates;
    set _gtf_label_candidates;
    if findw(upcase("&label_snps"),upcase(strip(SNP)),' ')>0;
  run;
  proc sort data=_gtf_label_candidates;
    by BP SNP;
  run;

  data _null_;
    set _gtf_label_candidates end=_last;
    length _ch $1 _positions $32767 _layout $12;
    array _labels[200] $128 _temporary_;
    array _bp[200] _temporary_;
    array _x[200] _temporary_;
    array _text_width[200] _temporary_;
    array _half[200] _temporary_;
    array _separation[200] _temporary_;
    array _centers[200] _temporary_;
    array _block_sum[200] _temporary_;
    array _block_count[200] _temporary_;
    array _block_first[200] _temporary_;
    array _block_last[200] _temporary_;
    retain _nlabels 0;
    _nlabels+1;
    if _nlabels>200 then do;
      put 'ERROR: SAS local-GTF label planner supports at most 200 targets.';
      abort cancel;
    end;
    _labels[_nlabels]=strip(SNP);
    _bp[_nlabels]=BP;
    _font=&font_size;
    _em=0;
    do _j=1 to lengthn(strip(SNP));
      _ch=substr(strip(SNP),_j,1);
      if indexc('WM@%',_ch)>0 then _em+0.95;
      else if indexc('ilI1.,:;|',_ch)>0 then _em+0.40;
      else if 'A'<=_ch and _ch<='Z' then _em+0.72;
      else _em+0.62;
    end;
    _text_width[_nlabels]=_em*_font*96/72*1.12+4;
    if not _last then return;

    _left=0.11*&design_width;
    _right=0.89*&design_width;
    _plot_width=_right-_left;
    _xmin=max(1,&center_bp-&window_bp);
    _xmax=&center_bp+&window_bp;
    _span=_xmax-_xmin;
    _expected=countw("&label_snps",' ');
    if _nlabels ne _expected or _span<=0 then do;
      put 'WARNING: SAS local-GTF planner lacks one or more target positions; using vertical SAS spacing.';
      call symputx("&out_angle",90,'g');
      call symputx("&out_positions",'','g');
      stop;
    end;
    do _i=1 to _nlabels;
      _x[_i]=_left+(_bp[_i]-_xmin)/_span*_plot_width;
    end;
    _request=lowcase("&requested_layout");
    _shift_limit=min(90,0.12*_plot_width);
    _chosen=0;
    do _mode=1 to 2 while(_chosen=0);
      if _mode=1 then _gap=_font*96/72*0.7+6;
      else _gap=5;
      _needed=(_nlabels-1)*_gap;
      do _i=1 to _nlabels;
        if _mode=1 then _half[_i]=_text_width[_i]/2;
        else _half[_i]=_font*96/72*0.72;
        _needed+2*_half[_i];
      end;
      if _needed>_plot_width then do;
        if _mode=1 and _request='horizontal' then do;
          put 'ERROR: Horizontal target labels cannot fit; use auto/vertical or increase gtf_design_width.';
          abort cancel;
        end;
        if _mode=2 then do;
          put 'ERROR: Too many vertical target labels; increase gtf_design_width.';
          abort cancel;
        end;
        continue;
      end;
      _separation[1]=0;
      do _i=2 to _nlabels;
        _separation[_i]=_separation[_i-1]+_half[_i-1]+_half[_i]+_gap;
      end;
      _blocks=0;
      do _i=1 to _nlabels;
        _blocks+1;
        _block_sum[_blocks]=_x[_i]-_separation[_i];
        _block_count[_blocks]=1;
        _block_first[_blocks]=_i;
        _block_last[_blocks]=_i;
        do while(_blocks>=2);
          if _block_sum[_blocks-1]/_block_count[_blocks-1]
            <=_block_sum[_blocks]/_block_count[_blocks] then leave;
          _block_last[_blocks-1]=_block_last[_blocks];
          _block_sum[_blocks-1]=_block_sum[_blocks-1]+_block_sum[_blocks];
          _block_count[_blocks-1]=_block_count[_blocks-1]+_block_count[_blocks];
          _blocks=_blocks-1;
        end;
      end;
      _minimum=_left+_half[1];
      _maximum=_right-_half[_nlabels]-_separation[_nlabels];
      do _b=1 to _blocks;
        _base=max(_minimum,min(_maximum,_block_sum[_b]/_block_count[_b]));
        do _i=_block_first[_b] to _block_last[_b];
          _centers[_i]=_base+_separation[_i];
        end;
      end;
      _max_shift=0;
      do _i=1 to _nlabels;
        _max_shift=max(_max_shift,abs(_centers[_i]-_x[_i]));
      end;
      if _mode=1 and _request='vertical' then continue;
      if _mode=1 and _request='auto' and _max_shift>_shift_limit then continue;
      _chosen=_mode;
    end;
    _positions='';
    do _i=1 to _nlabels;
      _adjusted_bp=round(_xmin+(_centers[_i]-_left)/_plot_width*_span,1);
      _positions=catx(' ',_positions,cats(strip(_labels[_i]),'=',put(_adjusted_bp,best32.)));
    end;
    _angle=ifn(_chosen=1,0,90);
    if _chosen=1 then _text_height=_font*96/72*1.2;
    else do;
      _text_height=0;
      do _i=1 to _nlabels;
        _text_height=max(_text_height,_text_width[_i]);
      end;
    end;
    _headroom_px=_text_height+max(12,_font*96/72*0.8);
    _headroom=max(0.035,min(0.45,_headroom_px/(&design_height*0.75)));
    _center_offset=_headroom*&design_height*0.75/(2*_text_height);
    call symputx("&out_angle",_angle,'g');
    call symputx("&out_positions",_positions,'g');
    call symputx("&out_headroom",_headroom,'g');
    call symputx("&out_center_offset",_center_offset,'g');
    put 'NOTE: SAS local-GTF planner chose angle=' _angle
        ' with ' _nlabels ' labels, headroom=' _headroom
        ' center_offset=' _center_offset ' and adjusted positions ' _positions;
  run;
%mend;
