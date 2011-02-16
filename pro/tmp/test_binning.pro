PRO test_binning_binbynum, data, w1, w2, w3, w4, ind1, ind2, equal=equal, minx=minx, maxx=maxx

  ;; currently only works for 4 bins
  IF keyword_set(equal) THEN BEGIN 
      perc1 = 0.25
      perc2 = 0.25
      perc3 = 0.25
      perc4 = 0.25
  ENDIF ELSE BEGIN 
      perc1 = 15250./21991.
      perc2 = 3840./21991.
      perc3 = 1940./21991.
      perc4 = 960./21991.
  ENDELSE 

  ndata = n_elements(data)
  w=lindgen(ndata)
  IF n_elements(minx) NE 0 THEN BEGIN 
      w2=where( data[w] GT minx, ndata)
      w=w[w2]
  ENDIF 
  IF n_elements(maxx) NE 0 THEN BEGIN 
      w2=where( data[w] LT maxx, ndata)
      w=w[w2]
  ENDIF 

  s=sort(data[w])
  w=w[s]

  bin1 = long( perc1*ndata )
  bin2 = long( perc2*ndata )
  bin3 = long( perc3*ndata )
  bin4 = long( perc4*ndata )

  w1 = w[0:bin1-1]
  w2 = w[bin1:bin1+bin2-1]
  w3 = w[bin1+bin2:bin1+bin2+bin3-1]
  w4 = w[bin1+bin2+bin3:ndata-1]

  ind1 = [0L,     bin1,        bin1+bin2         bin1+bin2+bin3]
  ind2 = [bin1-1, bin1+bin2-1, bin1+bin2+bin3-1, nw-1          ]


END 

PRO test_binning, yynorm, yypow, norm, pow, nodisplay=nodisplay, equal=equal,$
                  addnoise=addnoise

  ;; must have 2 windows started

  minx = 10.
  maxx = 90.

  x0 = 0.01
  x1 = 100
  
  nx = 1000L
  
;  xx_init = arrscl( randomu(seed, nx), x0, x1, arrmin=0., arrmax=1.)
;  xx_init = xx_init[sort(xx_init)]
  xx_init = arrscl( findgen(nx), x0, x1 )
  yy_init = yynorm*(xx_init)^yypow

  ;; vals for fit
  np=400
  nn=400

  IF keyword_set(addnoise) THEN BEGIN 
      ;;  xerrfac = 100.
      yerrfac = 0.1
      ;;  xerr = sqrt(xx_init > 0.)/xerrfac
      yerr = sqrt(yy_init > 0.)/yerrfac

      ;;xx = xx_init+xerr*randomn(seed, nx)
      yy = yy_init+yerr*randomn(seed, nx)
  ENDIF ELSE BEGIN 
      xx=xx_init
      yy=yy_init
      xerr=fltarr(nx)
      yerr=fltarr(nx)
  ENDELSE 

  IF NOT keyword_set(nodisplay) THEN BEGIN 
      ;wset,0
      ;xrange=1.1*prange(xx,xerr)
      ;yrange=1.1*prange(yy,yerr)
      ww=where(xx_init GT minx)
      ;xrange[0] = xrange[0] > 0.2*minx
      ;yrange[0] = yrange[0] > 0.2*yy_init[ww[0]]
      plot,xx,yy,psym=1;,xrange=xrange,yrange=yrange,$
        ;/xlog,/ylog
      oplot,xx_init,yy_init,color=!red
      rrr = 1.e8
      oplot, [minx,minx],[0.1,rrr],color=!yellow
      oplot, [maxx,maxx],[0.1,rrr],color=!yellow
  ENDIF 

  IF NOT keyword_set(addnoise) THEN BEGIN 
      xerr=replicate(1.0,nx)
      yerr=xerr
  ENDIF 
  test_binning_binbynum, xx, w1, w2, w3, w4, minx=minx, maxx=maxx, equal=equal

  ;inputweights = yy[w1]
  wmom, xx[w1], yerr[w1], meanx1, tsig, mxerr1, inputweights=inputweights,/noweight
  wmom, yy[w1], yerr[w1], meany1, tsig, myerr1, inputweights=inputweights,/noweight

  ;inputweights = yy[w2]
  wmom, xx[w2], yerr[w2], meanx2, tsig, mxerr2, inputweights=inputweights,/noweight
  wmom, yy[w2], yerr[w2], meany2, tsig, myerr2, inputweights=inputweights,/noweight

  ;inputweights = yy[w3]
  wmom, xx[w3], yerr[w3], meanx3, tsig, mxerr3, inputweights=inputweights,/noweight
  wmom, yy[w3], yerr[w3], meany3, tsig, myerr3, inputweights=inputweights,/noweight

  ;inputweights = yy[w4]
  wmom, xx[w4], yerr[w4], meanx4, tsig, mxerr4, inputweights=inputweights,/noweight
  wmom, yy[w4], yerr[w4], meany4, tsig, myerr4, inputweights=inputweights,/noweight

  mx =  [meanx1, meanx2, meanx3, meanx4]
  my =  [meany1, meany2, meany3, meany4]
  mxerr = [mxerr1, mxerr2, mxerr3, mxerr4]
  myerr = [myerr1, myerr2, myerr3, myerr4]
  
  IF NOT keyword_set(nodisplay) THEN BEGIN 
      mclr=!magenta
      oploterror, mx, my, mxerr, myerr, psym=4, color=mclr,errcolor=mclr, $
        thick=2.*!p.thick
  ENDIF 
  ;; do some fits
  xold = !x
  yold = !y
  pold = !p

  ;wset,1
  aguess = [yynorm, yypow]
  fitpower, mx, my, myerr, aguess, yfit, aout, aerr
  nsig = 200.

  nvhigh = ( aout[0]+nsig*aerr[0] ) < max([4.*yynorm,4.*aout[0]])
  IF nvhigh LT yynorm THEN nvhigh = yynorm*1.1
  nvlow = ( aout[0]- nsig*aerr[0] ) > min([0.2*yynorm,0.2*aout[0]])
  IF nvlow GT yynorm THEN nvlow = yynorm*0.9

  nphigh = ( aout[1]+nsig*aerr[1] ) < max([4.*yypow,4.*aout[1]])
  IF nphigh LT yypow THEN nphigh = yypow*1.1
  nplow = ( aout[1]- nsig*aerr[1] ) > min([0.2*yypow,0.2*aout[1]])
  IF nplow GT yypow THEN nplow = yypow*0.9

  normvals = arrscl( findgen(nn), nvlow, nvhigh )
  powvals = arrscl( findgen(np), nplow, nphigh )

  px1 = 0.37
  py1 = 0.7
  px2 = 0.65
  py2 = 0.93

  !p.position = [ [px1, py1], [px2, py2] ]
  pow_chisq_conf, mx, my, myerr, powvals, normvals, $
    chisq_surf, pow, norm, powlow, powhigh, normlow, normhigh,$
    powallow=powallow, normallow=normallow, nodisplay=nodisplay,/noerase,$
    ytitle='Norm',xtitle=!tsym.alpha
  IF NOT keyword_set(nodisplay) THEN $
    oplot,[yypow], [yynorm], psym=7, symsize=2, color=!blue
  range2error, powlow[0], pow, powhigh[1], peh, pel
  range2error, normlow[0], norm, normhigh[1], neh, nel

  colprint,['Norm  ','Power  '],ntostr(aout),' '+ntostr([norm,pow])+' +/- '+ntostr([max([neh,nel]),max([peh,pel])])
  print
;  ploterror, mx, my, mxerr, myerr, psym=1, /xlog, /ylog, $
;    xrange=[10, 100], yrange=[10, 100]
;  ploterror, mx, my, mxerr, myerr, psym=1

  ;;forprint,mx, my,myerr
;  oplot,xx_init,yy_init,color=!red

  s=sort(xx_init)
  yff = norm*xx_init[s]^pow
  ;wset,0
  !x=xold
  !y=yold
  !p=pold

  IF NOT keyword_set(nodisplay) THEN BEGIN 
      oplot, xx_init[s], yff, color=!blue
  ENDIF 

END 
