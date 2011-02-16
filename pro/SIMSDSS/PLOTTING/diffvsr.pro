PRO diffvsr, cat, drplot=drplot, plt2=plt2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;   Make plot of e1meas - e2real  vs r
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_params() EQ 0 THEN BEGIN
    print,'-Syntax: diffvsr, cat, wg, e1, e2'
    return
  ENDIF 

  tflags = cat[ rem_dup(cat.typeflag) ].typeflag
  print,'Found ',strtrim(string(n_elements(tflags)),2),'  distinct PSF shapes'

  FOR i=0, n_elements(tflags)-1 DO BEGIN
    w=where(cat.typeflag EQ tflags[i])
    
    IF keyword_set(drplot) THEN BEGIN
      correct, cat[w], e1, e2, err, wg=wg, psfe1=p1,psfe2=p2,/silent,/drplot
      key=get_kbrd(20)
      IF key EQ 'q' THEN return
    ENDIF ELSE BEGIN
      correct, cat[w], e1, e2, err, wg=wg, psfe1=p1,psfe2=p2,/silent
    ENDELSE 

    fitlin, cat[w[wg]].sr, e1[wg]-cat[w[wg]].se1, err[wg],$
      aa1,sigaa1,bb1,sigbb1,/silent
    fitlin, cat[w[wg]].sr, e2[wg]-cat[w[wg]].se2, err[wg],$
      aa2,sigaa2,bb2,sigbb2,/silent
    

    IF i EQ 0 THEN BEGIN
      a1=aa1
      siga1=sigaa1
      b1=bb1
      sigb1=sigbb1

      a2=aa2
      siga2=sigaa2
      b2=bb2
      sigb2=sigbb2

      psfe1=p1
      psfe2=p2
    ENDIF ELSE BEGIN 
      a1=[a1,aa1]
      siga1=[siga1,sigaa1]
      b1=[b1,bb1]
      sigb1=[sigb1,sigbb1]

      a2=[a2,aa2]
      siga2=[siga2,sigaa2]
      b2=[b2,bb2]
      sigb2=[sigb2,sigbb2]

      psfe1=[psfe1,p1]
      psfe2=[psfe2,p2]
    ENDELSE
  ENDFOR 

  pold=!p.multi
  !p.multi=[0,2,2]
  yrange=[-1,1]

  xtitle='e1 Psf'
  title='Offset of (e1 - e1 input) vs R line'
  ytitle='Offset'
  plot,psfe1,a1,psym=7,yrange=yrange,xtitle=xtitle,ytitle=ytitle,title=title
  oploterr,psfe1,a1,siga1
  oplot,[-10,10],[0,0]

  title='Slope of (e1 - e1 input) vs R line'
  ytitle='Slope'
  plot,psfe1,b1,psym=7,yrange=yrange,xtitle=xtitle,ytitle=ytitle,title=title
  oploterr,psfe1,b1,sigb1
  oplot,[-10,10],[0,0]

  xtitle='e1 Psf'
  title='Offset of (e2 - e2 input) vs R line'
  ytitle='Offset'
  plot,psfe1,a2,psym=7,yrange=yrange,xtitle=xtitle,ytitle=ytitle,title=title
  oploterr,psfe1,a2,siga2
  oplot,[-10,10],[0,0]

  title='Slope of (e2 - e2 input) vs R line'
  ytitle='Slope'
  plot,psfe1,b2,psym=7,yrange=yrange,xtitle=xtitle,ytitle=ytitle,title=title
  oploterr,psfe1,b2,sigb2
  oplot,[-10,10],[0,0]

  key=get_kbrd(20)
  IF key EQ 'q' THEN BEGIN
    return & !p.multi=pold
  ENDIF 

  xtitle='e2 Psf'
  title='Offset of (e1 - e1 input) vs R line'
  ytitle='Offset'
  plot,psfe2,a1,psym=7,yrange=yrange,xtitle=xtitle,ytitle=ytitle,title=title
  oploterr,psfe1,a1,siga1
  oplot,[-10,10],[0,0]

  title='Slope of (e1 - e1 input) vs R line'
  ytitle='Slope'
  plot,psfe2,b1,psym=7,yrange=yrange,xtitle=xtitle,ytitle=ytitle,title=title
  oploterr,psfe1,b1,sigb1
  oplot,[-10,10],[0,0]

  xtitle='e2 Psf'
  title='Offset of (e2 - e2 input) vs R line'
  ytitle='Offset'
  plot,psfe2,a2,psym=7,yrange=yrange,xtitle=xtitle,ytitle=ytitle,title=title
  oploterr,psfe1,a2,siga2
  oplot,[-10,10],[0,0]

  title='Slope of (e2 - e2 input) vs R line'
  ytitle='Slope'
  plot,psfe2,b2,psym=7,yrange=yrange,xtitle=xtitle,ytitle=ytitle,title=title
  oploterr,psfe1,b2,sigb2
  oplot,[-10,10],[0,0]


  !p.multi=pold

return
END







