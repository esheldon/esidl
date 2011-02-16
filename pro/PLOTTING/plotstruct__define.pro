function plotstruct::init
  return,1
end 

pro correlate::pmulti, struct, xtagname, ytagname, yerrtagname, $
             nsplit=nsplit, $
             xtitle=xtitle, ytitle=ytitle, $
             xrange=xrange, yrange=yrange, $
             xstyle=xstyle, ystyle=ystyle, $
             _extra=_extra

  if not tag_exist(struct, xtagname, index=xtag) then message,'X Tag does not exist: '+ntostr(xtagname)
  if not tag_exist(struct, ytagname, index=ytag) then message,'Y Tag does not exist: '+ntostr(ytagname)
  if not tag_exist(struct, yerrtagname, index=yerrtag) then message,'Y Err Tag does not exist: '+ntostr(yerrtagname)

  if n_elements(nsplit) ne 0 then begin
      ;; this means that we should overplot the splits.
      ;; the actual number of plots will be nbin/nsplit
      nbin = n_elements(struct)/nsplit
      psyms = -[8, 4, 1, 5, 6, 7]
      npsym = n_elements(psyms)
      linestyles=[0,1,3,5]
      nline = n_elements(linestyles)
  endif else begin 
      nsplit = 1
      nbin = n_elements(struct)
      psym = 8
  endelse 

  simpctable, colorlist=colorlist

  nc = n_elements(colorlist)

  if nbin gt 1 then begin 
      mplot_value = self->mplot_value(nbin)
      erase & multiplot, mplot_value, /square
  endif else begin 
      aspect = 1
  endelse 


  lab = self->axis_labeled(nbin)
  plotstruct = self->plotstruct(sharr, type, $
                                xlog=xlog, xmin=xmin, xmax=xmax, $
                                ylog=ylog, ymin=ymin, ymax=ymax, $
                                symmetric=symmetric)


  if nbin gt 8 and keyword_set(dops) then begin 
      !x.thick = 2
      !y.thick = 2
      !p.thick = 2
  endif 

  j = 0
  ish = 0
  for i=0l, nbin-1 do begin 

      for j=0, nsplit-1 do begin 

          if j gt 0 then overplot=1 else overplot=0

          if keyword_set(incolor) then begin 
              color=colorlist[j mod nc]
              add_arrval, color, lcolors
          endif 
          if n_elements(psyms) ne 0 then begin 
              psym = psyms[j mod npsym]
              lpsym = psym
              if nsplit ne 1 then lpsym=abs(lpsym)
              add_arrval, lpsym, lpsyms
          endif 
          if n_elements(linestyles) ne 0 then begin 
              linestyle=linestyles[j mod nline]
              add_arrval, linestyle, llinestyles
          endif 
          self->_plot_profile, sharr[ish], plotstruct, $
            label=label, $
            xlabeled=lab.xlabeled[i], $
            ylabeled=lab.ylabeled[i], $
            overplot=overplot, $
            aspect=aspect, psym=psym, linestyle=linestyle, $
            color=color

          if n_elements(labels) ne 0 then begin 
              add_arrval, labels[ish], llabels
          endif 

          ish = ish+1
      endfor 

      if n_elements(llabels) ne 0 then begin 
          ;;help,llabels,llinestyles,lpsyms,lcolors
          if n_elements(label_charsize) eq 0 then label_charsize=1

          legend, llabels, line=llinestyles, psym=lpsyms, color=lcolors, $
            /right, box=0, charsize=label_charsize, margin=0
      endif 
      
      delvarx, lcolors, lpsyms, llinestyles, llabels

      if nbin gt 1 and i ne nbin-1 then multiplot
  endfor 


  if nbin gt 1 then begin 
      multiplot, /reset
  endif 

  if n_elements(charsize) ne 0 then !p.charsize=chold

  if keyword_set(dops) then begin 
      if keyword_set(encapsulated) then trim_bbox=1
      if keyword_set(landscape) then landfix=1
      endplot,trim_bbox=trim_bbox, landfix=landfix
  endif 


end 


FUNCTION plotstruct::plotstruct, sh, type, $
                   xlog=xlog, xmin=xmin, xmax=xmax, $
                   ylog=ylog, ymin=ymin, ymax=ymax, $
                   symmetric=symmetric
  
  nbin = n_elements(sh)

  ;; Tag info
  tg=self->typetags(sh, type)
  tit = self->titles(type)

  data = sh.(tg.tag)
  err = sh.(tg.errtag)
  IF type EQ 'clustcorr' THEN BEGIN 
      data = data - 1.0
  ENDIF 

  yrange = prange( data, err, slack=1.5, symmetric=symmetric)
  xrange = [min(sh.meanr), max(sh.meanr)]/1000

  IF keyword_set(xlog) THEN BEGIN 
      xrange[0] = 10.0/1000.0   ; 10 kpc
 

      xrange[1] = xrange[1]*2.0
      xtickformat = 'loglabels'
      xticklen = 0.04
  ENDIF ELSE BEGIN 
      xlog=0
      xtickformat = ''
      xticklen = !x.ticklen
  ENDELSE 
  IF n_elements(xmin) NE 0 THEN xrange[0] = xmin
  IF n_elements(xmax) NE 0 THEN xrange[1] = xmax

  IF keyword_set(ylog) THEN BEGIN 
      yrange[0] = 1.e-2
      yrange[1] = yrange[1]*1.5
      ytickformat = 'loglabels'
      yticklen = 0.04
  ENDIF ELSE BEGIN 
      ylog=0
      ytickformat = ''
      yticklen = !y.ticklen
  ENDELSE 
  IF n_elements(ymin) NE 0 THEN yrange[0] = ymin
  IF n_elements(ymax) NE 0 THEN yrange[1] = ymax

  st = $
    {type: type, $
     tag: tg.tag, $
     errtag: tg.errtag, $
     $
     xlog: xlog, $
     xtickformat: xtickformat, $
     xticklen: xticklen, $
     xrange: xrange, $
     xtitle: tit.xtitle, $
     $
     ylog: ylog, $
     ytickformat: ytickformat, $
     yticklen: yticklen, $
     yrange: yrange, $
     ytitle: tit.ytitle $
    }  

  return,st
END 



function plotstruct::mplot_value, nbin

  case nbin of
      1: message,"nbin=1 doesn't make sense for multiplot"
      2: mplot_value=[1,2]
      3: mplot_value=[3,1]
      6: mplot_value=[3,2]
      8: mplot_value=[4,2]
      9: mplot_value=[3,3]
      12: mplot_value=[4,3]
      16: mplot_value=[4,4]
      else: message,'unknown nbin: '+ntostr(nbin)
  endcase 
  return,mplot_value

end 

function plotstruct::axis_labeled, nbin
  case nbin of
      1: begin 
          xlabeled = 1
          ylabeled = 1
      end 
      2: begin 
          xlabeled = [0,1]
          ylabeled = [1,1]
      end 
      3: begin 
          xlabeled = [1,1,1]
          ylabeled = [1,0,0]
      end 
      8: begin 
          xlabeled = [0,0,0,0,1,1,1,1]
          ylabeled = [1,0,0,0,1,0,0,0]
      end 
      12: begin 
          xlabeled = [0,0,0,0,0,0,0,0,1,1,1,1]
          ylabeled = [1,0,0,0,1,0,0,0,1,0,0,0]
      end 
      16: begin 
          xlabeled = [0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1]
          ylabeled = [1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0]
      end 
      else: message,'unsupported nbin: '+ntostr(nbin)
  endcase 
  st = {xlabeled: xlabeled, $
        ylabeled: ylabeled}
  return,st
end 


pro plotstruct__define
  struct = { $
             plotstruct, $
             plostruct_dummyvar: 0 $
           }
end 
