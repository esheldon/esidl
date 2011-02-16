PRO corrplt, cat, norho=norho, psfile=psfile, median=median

  IF n_params() EQ 0 THEN BEGIN
    print,'-Syntax: corrplt, cat, norho=norho, psfile=psfile, median=median'
    return
  ENDIF 
 
  IF keyword_set(norho) THEN norho = 1 ELSE norho = 0
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Find out how many different psf's we have.  This is determined by
  ; how many different typeflags there are in cat.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  tflags = cat[ rem_dup(cat.typeflag) ].typeflag
  print,'Found ',strtrim(string(n_elements(tflags)),2),'  distinct PSF shapes'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Loop over the different psf's
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ngal = 0
  FOR i=0, n_elements(tflags)-1 DO BEGIN
    w = where(cat.typeflag EQ tflags[i])
    correct, cat[w], e1, e2, err, wg=wg, psfe1=p1,psfe2=p2,$
             e1mom=e1mom,e2mom=e2mom, olde1mom=olde1mom, olde2mom=olde2mom, $
             norho=norho, /silent
    
    ;; first element is weighted mean, second weighted sigma, third is median
    IF keyword_set(median) THEN BEGIN
      m1 = e1mom[2]
      m2 = e2mom[2]
      om1 = olde1mom[2]
      om2 = olde2mom[2]
    ENDIF ELSE BEGIN
      m1 = e1mom[0]
      m2 = e2mom[0]
      om1 = olde1mom[0]
      om2 = olde2mom[0]
    ENDELSE 
    IF i EQ 0 THEN BEGIN
      me1 = m1
      me2 = m2             ;;;  corrected stuff
      sige1 = e1mom[1]
      sige2 = e2mom[1]

      ome1 = om1
      ome2 = om2           ;;; uncorrected stuff
      osige1 = olde1mom[1]
      osige2 = olde2mom[1]

      psfe1 = p1
      psfe2 = p2
    ENDIF ELSE BEGIN
      me1 = [ me1,m1 ]
      me2 = [ me2,m2 ]           ;;;  corrected stuff
      sige1 = [ sige1,e1mom[1] ]
      sige2 = [ sige2,e2mom[1] ]

      ome1 = [ome1, om1]
      ome2 = [ome2, om2]         ;;; uncorrected stuff
      osige1 = [ osige1, olde1mom[1] ]
      osige2 = [ osige2, olde2mom[1] ]

      psfe1 = [ psfe1,p1 ]
      psfe2 = [ psfe2,p2 ]
    ENDELSE    
    ngal = ngal + n_elements(wg)
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; do some plots
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  myusersym, 'fill_circle'
  psym=8
  
  IF NOT keyword_set(yrange) THEN yrange=[-1,1]
  xrange=[-.5,.5]
  xvals = findgen(100)/10.0-5.0  ;;plenty big

  pold=!p.charsize
  !p.charsize=1.0

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Uncorrected stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  yrange=prange(ome1,osige1,ome2,osige2,/slack)
  
  erase & multiplot, [2,2], /square
  
  fitlin, psfe1, ome1, osige1, a, siga, b, sigb,/silent
  out = 'a = '+strmid(strtrim(string(a),2),0,6)+'  b = ' $
    +strmid(strtrim(string(b),2),0,6) 
  ploterror, psfe1, ome1, osige1, ytitle='e!D1!N Gal',$
             psym=psym,yrange=yrange,xrange=xrange,ystyle=1,xstyle=1,$
             charsize=1.0,title=strtrim(string(ngal),2)+'  Uncorrected Galaxies'
  oplot,[-10,10],[0,0]
  oplot,xvals,a + b*xvals, color=!blue
  xyouts,.9*xrange[0],.85*yrange[1],out
  
  multiplot
  fitlin, psfe2, ome1, osige1, a, siga, b, sigb,/silent
  out = 'a = '+strmid(strtrim(string(a),2),0,6)+'  b = ' $
    +strmid(strtrim(string(b),2),0,6)
  ploterror, psfe2, ome1, osige1,$
             psym=psym,yrange=yrange,xrange=xrange,ystyle=1,xstyle=1
  oplot,[-10,10],[0,0]
  oplot,xvals,a + b*xvals, color=!blue
  ;;xyouts,.9*xrange[0],.85*yrange[1],out
  
  multiplot
  fitlin, psfe1, ome2, osige2, a, siga, b, sigb,/silent
  out = 'a = '+strmid(strtrim(string(a),2),0,6)+'  b = ' $
    +strmid(strtrim(string(b),2),0,6)
  ploterror, psfe1, ome2, osige2, xtitle='e!D1!N PSF', ytitle='e!D2!N Gal',$
             psym=psym,yrange=yrange,xrange=xrange,ystyle=1,xstyle=1
  oplot,[-10,10],[0,0]
  oplot,xvals,a + b*xvals, color=!blue
  ;;xyouts,.9*xrange[0],.85*yrange[1],out
  
  multiplot
  fitlin, psfe2, ome2, osige2, a, siga, b, sigb,/silent
  out = 'a = '+strmid(strtrim(string(a),2),0,6)+'  b = ' $
    +strmid(strtrim(string(b),2),0,6)
  ploterror, psfe2, ome2, xtitle='e!D2!N PSF',$
             psym=psym,yrange=yrange,osige2, xrange=xrange,ystyle=1,xstyle=1
  oplot,[-10,10],[0,0]
  oplot,xvals,a + b*xvals, color=!blue
  xyouts,.9*xrange[0],.85*yrange[1],out
  
;  x=!d.x_size/2.0
;  y=.97*!d.y_size
;  xyouts,x,y,/device,align=.5,strtrim(string(ngal),2)+'  Uncorrected Galaxies'
  
  multiplot,/reset
  
  IF display_type() EQ 'X' THEN key = get_kbrd(20)
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Corrected stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF norho THEN tt='First order correction' $
  ELSE tt = 'Higher order correction'

  erase & multiplot, [2,2], /square
  
  fitlin, psfe1, me1, sige1, a, siga, b, sigb,/silent
  out = 'a = '+strmid(strtrim(string(a),2),0,6)+'  b = ' $
    +strmid(strtrim(string(b),2),0,6)
  ploterror, psfe1, me1, sige1, ytitle='e!D1!N Gal',$
             psym=psym,yrange=yrange,xrange=xrange,ystyle=1,xstyle=1,$
             title=strtrim(string(ngal),2)+'  Corrected Galaxies'
  oplot,[-10,10],[0,0]
  oplot,xvals,a + b*xvals, color=!blue
  xyouts,.9*xrange[0],.85*yrange[1],out
  
  multiplot
  fitlin, psfe2, me1, sige1, a, siga, b, sigb,/silent
  out = 'a = '+strmid(strtrim(string(a),2),0,6)+'  b = ' $
    +strmid(strtrim(string(b),2),0,6)
  ploterror, psfe2, me1, sige1,$
             psym=psym,yrange=yrange,xrange=xrange,ystyle=1,xstyle=1,$
             title=tt
  oplot,[-10,10],[0,0]
  oplot,xvals,a + b*xvals, color=!blue
  xyouts,.9*xrange[0],.85*yrange[1],out
  
  multiplot
  fitlin, psfe1, me2, sige2, a, siga, b, sigb,/silent
  out = 'a = '+strmid(strtrim(string(a),2),0,6)+'  b = ' $
    +strmid(strtrim(string(b),2),0,6)
  ploterror, psfe1, me2, sige2, xtitle='e!D1!N PSF', ytitle='e!D2!N Gal',$
             psym=psym,yrange=yrange,xrange=xrange,ystyle=1,xstyle=1
  oplot,[-10,10],[0,0]
  oplot,xvals,a + b*xvals, color=!blue
  xyouts,.9*xrange[0],.85*yrange[1],out
  
  multiplot
  fitlin, psfe2, me2, sige2, a, siga, b, sigb,/silent
  out = 'a = '+strmid(strtrim(string(a),2),0,6)+'  b = ' $
    +strmid(strtrim(string(b),2),0,6)
  ploterror, psfe2, me2, sige2, xtitle='e!D2!N PSF',$
             psym=psym,yrange=yrange,xrange=xrange,ystyle=1,xstyle=1
  oplot,[-10,10],[0,0]
  oplot,xvals,a + b*xvals, color=!blue
  xyouts,.9*xrange[0],.85*yrange[1],out
  
  
  multiplot,/reset

  !p.charsize=pold

return 
END






