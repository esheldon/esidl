PRO corshape_again_getrange, x1, y1, x2, y2, nperbin, xrange, yrange

  binner_bynum, x1, y1, nperbin, txo1, tyo1, tsig1, tnum1
  binner_bynum, x2, y2, nperbin, txo2, tyo2, tsig2, tnum2

  w1=where(tnum1 EQ nperbin)
  w2=where(tnum2 EQ nperbin)
  xrange=prange(txo1[w1],txo2[w2],/noerror,slack=0.3,/symmetric)
  yrange=prange(tyo1[w1],tyo2[w2],/noerror,slack=0.4,/symmetric)

END 

PRO corrforslope, ac, newe1gal, newe2gal

  newe1gal = ac.gale1_corrected - $
    (ac.intercept_ge1pe1 + ac.slope_ge1pe1*ac.galpsfe1)
  newe2gal = ac.gale2_corrected - $
    (ac.intercept_ge2pe2 + ac.slope_ge2pe2*ac.galpsfe2)

  ;newe1gal = ac.gale1_corrected - $
  ;  (ac.intercept_ge1pe1R + ac.slope_ge1pe1R*ac.galpsfe1)
  ;newe2gal = ac.gale2_corrected - $
  ;  (ac.intercept_ge2pe2R + ac.slope_ge2pe2R*ac.galpsfe2)

END 

PRO plot_egal_vs_epsf, ac, copy=copy, title=title

  ocolor=!red
  fcolor=!blue
  xx=[-1.,1.]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plots of egal vs epsf
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;; about 8-9 bins
  frac = 1./9.
  nperbin = long( round(frac*n_elements(ac.galpsfe1)) )

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find out xrange for gals
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  corshape_again_getrange, ac.galpsfe1, ac.gale1, ac.galpsfe2, ac.gale2, $
                           nperbin, xrange, yrange

  erase & multiplot,[2,2],/square

  plot_fitlin, ac.galpsfe1,ac.gale1_corrected, nperbin, $
    int_ge1pe1, interr_ge1pe1, slope_ge1pe1, slopeerr_ge1pe1,$
    ytitle='e!D1!N Gal',xrange=xrange,yrange=yrange,$
    ystyle=1,xstyle=1, title=title

  cyrange=!y.crange 
  plothist,ac.galpsfe1,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, ac.galpsfe2,ac.gale1_corrected, nperbin, $
    int_ge1pe2, interr_ge1pe2, slope_ge1pe2, slopeerr_ge1pe2,$
    ystyle=8+1,xstyle=1,xrange=xrange,yrange=yrange
  axis, yaxis=1, ystyle=1, ytitle='e!D1!N Gal'

  cyrange=!y.crange 
  plothist,ac.galpsfe2,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, ac.galpsfe1,ac.gale2_corrected,nperbin, $
    int_ge2pe1, interr_ge2pe1, slope_ge2pe1, slopeerr_ge2pe1,$
    ytitle='e!D2!N Gal',xtitle='e!D1!N PSF',xstyle=1,ystyle=1,$
    xrange=xrange,yrange=yrange

  multiplot
  plot_fitlin, ac.galpsfe2,ac.gale2_corrected,nperbin, $
    int_ge2pe2, interr_ge2pe2, slope_ge2pe2, slopeerr_ge2pe2,$
    xtitle='e!D2!N PSF',$
    ystyle=8+1,xstyle=1,xrange=xrange,yrange=yrange
  axis,yaxis=1,ystyle=1, ytitle='e!D2!N Gal'

  multiplot, /reset

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plots of egal vs smear polarizability*epsf
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  corshape_again_getrange, ac.galpsfe1*ac.corr, ac.gale1, $
                           ac.galpsfe2*ac.corr, ac.gale2, $
                           nperbin, xrange, yrange

  erase & multiplot,[2,2],/square

  plot_fitlin, ac.galpsfe1*ac.corr,ac.gale1_corrected, nperbin, $
    int_ge1pe1R, interr_ge1pe1R, slope_ge1pe1R, slopeerr_ge1pe1R,$
    ytitle='e!D1!N Gal',xrange=xrange,yrange=yrange,$
    ystyle=1,xstyle=1, title=title

  cyrange=!y.crange 
  plothist,ac.galpsfe1*ac.corr,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, ac.galpsfe2*ac.corr,ac.gale1_corrected, nperbin, $
    int_ge1pe2R, interr_ge1pe2R, slope_ge1pe2R, slopeerr_ge1pe2R,$
    ystyle=8+1,xstyle=1,xrange=xrange,yrange=yrange
  axis, yaxis=1, ystyle=1, ytitle='e!D1!N Gal'

  cyrange=!y.crange 
  plothist,ac.galpsfe2*ac.corr,xhist,yhist,bin=ebin,/noplot
  nyhist=yhist*abs(cyrange[0])/max(yhist)*0.5 - abs(cyrange[0])
  oplot, xhist, nyhist, psym=10, color=!grey50

  multiplot
  plot_fitlin, ac.galpsfe1*ac.corr,ac.gale2_corrected,nperbin, $
    int_ge2pe1R, interr_ge2pe1R, slope_ge2pe1R, slopeerr_ge2pe1R,$
    ytitle='e!D2!N Gal',xtitle='e!D1!N Psf*R',xstyle=1,ystyle=1,$
    xrange=xrange,yrange=yrange

  multiplot
  plot_fitlin, ac.galpsfe2*ac.corr,ac.gale2_corrected,nperbin, $
    int_ge2pe2R, interr_ge2pe2R, slope_ge2pe2R, slopeerr_ge2pe2R,$
    xtitle='e!D2!N Psf*R',$
    ystyle=8+1,xstyle=1,xrange=xrange,yrange=yrange
  axis,yaxis=1,ystyle=1, ytitle='e!D2!N Gal'

  multiplot, /reset

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; egal vs smear polarizability
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ocolor=!red
  fcolor=!blue
  xx=[-1.,1.]
  xrange=[0,1]
  yrange=[-0.05,0.05]

  erase & multiplot,[1,2],/square

  nperbin2=long(round(nperbin/2.))

  plot_fitlin, ac.corr, ac.gale1_corrected, nperbin2,$
    int_ge1R, interr_ge1R, slope_ge1R, slopeerr_ge1R,$
    ystyle=1,xstyle=1,xrange=xrange,yrange=yrange,$
    ytitle='e!D1!N Gal',title=title
  multiplot

  plot_fitlin, ac.corr, ac.gale2_corrected, nperbin2,$
    int_ge2R, interr_ge2R, slope_ge2R, slopeerr_ge2R,$
    ystyle=1,xstyle=1,xrange=xrange,yrange=yrange,$
    ytitle='e!D2!N Gal',xtitle='R'
  
  multiplot,/reset

  IF keyword_set(copy) THEN BEGIN 
      ;; copy e vs epsf
      ac.intercept_ge1pe1=int_ge1pe1 & ac.intercept_ge1pe1_error=interr_ge1pe1
      ac.slope_ge1pe1=slope_ge1pe1   & ac.slope_ge1pe1_error=slopeerr_ge1pe1

      ac.intercept_ge2pe1=int_ge2pe1 & ac.intercept_ge2pe1_error=interr_ge2pe1
      ac.slope_ge2pe1=slope_ge2pe1   & ac.slope_ge2pe1_error=slopeerr_ge2pe1

      ac.intercept_ge1pe2=int_ge1pe2 & ac.intercept_ge1pe2_error=interr_ge1pe2
      ac.slope_ge1pe2=slope_ge1pe2   & ac.slope_ge1pe2_error=slopeerr_ge1pe2

      ac.intercept_ge2pe2=int_ge2pe2 & ac.intercept_ge2pe2_error=interr_ge2pe2
      ac.slope_ge2pe2=slope_ge2pe2   & ac.slope_ge2pe2_error=slopeerr_ge2pe2

      ;;  copy e vs r*epsf
      ac.intercept_ge1pe1R=int_ge1pe1R & ac.intercept_ge1pe1_error=interr_ge1pe1
      ac.slope_ge1pe1R=slope_ge1pe1R   & ac.slope_ge1pe1_error=slopeerr_ge1pe1

      ac.intercept_ge2pe1R=int_ge2pe1R & ac.intercept_ge2pe1_error=interr_ge2pe1
      ac.slope_ge2pe1R=slope_ge2pe1R   & ac.slope_ge2pe1_error=slopeerr_ge2pe1

      ac.intercept_ge1pe2R=int_ge1pe2R & ac.intercept_ge1pe2_error=interr_ge1pe2
      ac.slope_ge1pe2R=slope_ge1pe2R   & ac.slope_ge1pe2_error=slopeerr_ge1pe2

      ac.intercept_ge2pe2R=int_ge2pe2R & ac.intercept_ge2pe2_error=interr_ge2pe2
      ac.slope_ge2pe2R=slope_ge2pe2R   & ac.slope_ge2pe2_error=slopeerr_ge2pe2

      ;; copy e vs 
      ac.intercept_ge1R = int_ge1R & ac.intercept_ge1R_error = interr_ge1R
      ac.slope_ge1R = slope_ge1R & ac.slope_ge1R_error = slopeerr_ge1R

      ac.intercept_ge2R = int_ge2R & ac.intercept_ge2R_error = interr_ge2R
      ac.slope_ge2R = int_ge2R & ac.slope_ge2R_error = slopeerr_ge2R

  ENDIF 


END 

PRO corshape_again, acfile

  ;basedir='/sdss3/data4/corr756/1/objcs/'
  ;IF n_elements(ac) EQ 0 THEN ac=mrdfits(basedir+'1/corshape_000756_1_i_N1.fit',1)

  ac_corr_file = repstr(acfile, 'corshape', 'corshapeagain')
  ac_corr_psfile = repstr(ac_corr_file, '.fit','.ps')

  print,ac_corr_file,ac_corr_psfile

  ac = mrdfits(acfile, 1)

  begplot,name=ac_corr_psfile, /color
  !p.charsize=1.0

  plot_egal_vs_epsf, ac, title='corrected'

  ;key=get_kbrd(1)

  corrforslope, ac, newe1gal, newe2gal
  ac_corr=ac
  ac_corr.gale1_corrected = newe1gal
  ac_corr.gale2_corrected = newe2gal
  plot_egal_vs_epsf, ac_corr, /copy, title='re-corrected'

  endplot

  mwrfits, ac_corr, ac_corr_file, /create

END 
