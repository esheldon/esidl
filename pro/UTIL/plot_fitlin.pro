PRO plot_fitlin, x, y, nperbin, intercept, intercepterr, slope, slopeerr, $
                 yin_err=yin_err, $
                 legend_charsize=legend_charsize, $
                 _extra=extra, title=title, plotnum=plotnum, aspect=aspect

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: plot_fitlin, x, y, nperbin, intercept, intercepterr, slope, slopeerr, _extra=extra, title=title, plotnum=plotnum, aspect=aspect'
      return
  ENDIF 

  myusersym, 'fill_circle'
  psym=8

  ocolor=!p.color

  IF !d.name EQ 'PS' THEN BEGIN 
      fcolor = !blue
  ENDIF ELSE BEGIN 
      fcolor = !green
  ENDELSE 

  maxx = max( abs(x) )
  xx=[-10.0*maxx,10.0*maxx]

  binner_bynum, x, y, nperbin, xo, yo, sig, num, yin_err=yin_err
  ;; only use bins with nperbin
  w=where(num EQ nperbin, nw)
  IF nw GT 1 THEN BEGIN 
      xo=xo[w]
      yo=yo[w]
      sig=sig[w]
  ENDIF ELSE BEGIN 
      print,'Not enough bins!'
  ENDELSE  

  IF n_elements(aspect) EQ 0 THEN BEGIN 
      ploterror, xo, yo, sig, _extra=extra, title=title, psym=psym
  ENDIF ELSE BEGIN 
      aploterror, aspect, xo, yo, sig, _extra=extra, title=title, psym=psym
  ENDELSE 
  fitlin, xo, yo, sig, intercept, intercepterr, slope, slopeerr, /silent
  oplot,[-1,1],[0,0],color=ocolor
  
  oplot,xx,intercept + slope*xx, color=fcolor

  sa=strmid(strtrim(string(intercept,format='(e8.1)'),2),0,8)
;  ssiga=strmid(strtrim(string(abs(intercepterr/intercept)),2),0,5)
  ssiga=strmid(strtrim(string(intercepterr,format='(e8.1)'),2),0,8)

  sb=strmid(strtrim(string(slope,format='(e8.1)'),2),0,8)
;  ssigb=strmid(strtrim(string(abs(slopeerr/slope)),2),0,5)
  ssigb=strmid(strtrim(string(slopeerr,format='(e8.1)'),2),0,8)

;  stta='a='+sa+' '+ssiga
;  sttb='b='+sb+' '+ssigb
  stta='a='+sa+!csym.plusminus+ssiga
  sttb='b='+sb+!csym.plusminus+ssigb

  ;;!p.background=!white
  legend, [stta,sttb], /right, box=0, charsize=legend_charsize

  IF keyword_set(plotnum) THEN BEGIN 
      yrange=!y.crange

      ;; want zero to be at bottom of plot
      pnum = num*abs(yrange[0])/max(num)*0.9 - abs(yrange[0]) 

      forprint,num,pnum
      oplot,xo,pnum,psym=10,color=!grey50
  ENDIF 

END 
