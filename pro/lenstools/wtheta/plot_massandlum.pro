PRO plot_massandlum, logplot=logplot


  ;; This plots the "density" used in calc_m2l_omega and the 
  ;; luminosity density scaled by M/L

  pold=!p.multi

  !p.multi=[0,1,2]

  mstr=mrdfits('/sdss5/data0/lensout/stripe10/main_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_comb_corr_N2.fit',1)

;  mstr=mrdfits('/sdss5/data0/lensout/stripe10/main_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_comb_corr_N1.fit',1)

;  mstr=mrdfits('/sdss5/data0/lensout/stripe10/main_zgal_gal_stripe10_comb_corr_N2.fit',1)

  lstr=mrdfits('/sdss6/data0/wtheta/wthetalumweq_stripe10_rw_N1.fit',1)
  rlstr=mrdfits('/sdss6/data0/wtheta/wthetarandlumweq_stripe10_rw_N1.fit',1)

  wl=where(lstr.meanlum GT 0.0 AND lstr.meanr GT 200.)
  wm=where(mstr.meanr GT 200.)

  meanlum = lstr.meanlum-rlstr.meanlum
  meanlumerr = sqrt( lstr.meanlumerr^2 + rlstr.meanlumerr^2)

  alpha = 0.773
  fac = (2.-alpha)/alpha
  errfac = 1.32
  m2l = 292.095

  xtitle='Projected Radius (h!U'+!tsym.minus+'1!N kpc)'
;  ytitle = '!S'+!tsym.sigma_cap+'!R!A'+!tsym.minus+'!N ('+!tsym.ltequal+'R) '+$
;    !tsym.minus+' !S'+!tsym.sigma_cap+'!R!A'+!tsym.minus+$
;    '!N (R) (h M'+sunsymbol()+' pc!U'+!tsym.minus+'2!N)'+$
;    '* (2.-'+!tsym.alpha+')/'+!tsym.alpha

  ytitle='Density (M'+sunsymbol()+' pc!U-2!N)'

  yrange=[0,30]
  xrange=[0,2000]
  IF keyword_set(logplot) THEN BEGIN 
      ylog=1
      xlog=1
      yrange=[1,30]
      xrange=[100,3000]
  ENDIF 

  ploterror,mstr.meanr[wm],mstr.sigma[wm]*fac,mstr.sigmaerr[wm]*fac*errfac,$
    psym=1,xrange=xrange,yrange=yrange, xtitle=xtitle,ytitle=ytitle,$
    xlog=xlog,ylog=ylog,xstyle=1,ystyle=1
  oploterror,lstr.meanr[wl],meanlum[wl]/100.*m2l,meanlumerr[wl]/100.*m2l,psym=4,color=!red,errcolor=!red

  legend,['Mass density','Luminosity Density * M/L'],psym=[1,4],$
    color=[!p.color,!red],thick=[!p.thick,!p.thick], box=0,/right

  !p.multi=pold

END 
