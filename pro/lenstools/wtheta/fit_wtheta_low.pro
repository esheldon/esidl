PRO fit_wtheta_low, meanr, diffdensecomb, errorcomb, cumul, xx, yycomb

  setupplot, dtype


  type='low_wtheta'

  stripestr = ['stripe10','stripe82','stripe36','stripe37','stripe42','stripe43','all']
  nstripe = n_elements(stripestr) - 1

  indir='/sdss5/data0/lensout/'+stripestr+'/'
  files = strarr(nstripe)
  rfiles = files
  FOR i=0L, nstripe-1 DO BEGIN 
      files[i] = indir[i]+type+'_'+stripestr[i]+'_sum_N1.fit'
      rfiles[i] = indir[i]+type+'rand_'+stripestr[i]+'_sum_N1.fit'
  ENDFOR 

  ;; first 2 (10,82) are looked at separately
  combine_wtheta,files[0],wt
  combine_wtheta,rfiles[0],wtr
  combine_wtheta,files[1],wt2
  combine_wtheta,rfiles[1],wtr2
  combine_wtheta,files,wtcomb
  combine_wtheta,rfiles,wtrcomb


  w=where(wt.npair NE 0.,nbin)
  w2=where(wt2.npair NE 0.,nbin2)
w=w2
  ;; density in #/Mpc^2
  diffdense = ( (wt.npair/wt.nlenses)-(wtr.npair/wtr.nlenses))/wt.area*1000.^2
  diffdense2 = ( (wt2.npair/wt2.nlenses)-(wtr2.npair/wtr2.nlenses))/wt2.area*1000.^2
  diffdensecomb = ( (wtcomb.npair/wtcomb.nlenses)-(wtrcomb.npair/wtrcomb.nlenses))/wtcomb.area*1000.^2

  error = sqrt( wt.npair/wt.nlenses^2 $
                + wtr.npair/wtr.nlenses^2 )/wt.area*1000.^2
  error2 = sqrt( wt2.npair/wt2.nlenses^2 $
                + wtr2.npair/wtr2.nlenses^2 )/wt2.area*1000.^2
  errorcomb = sqrt( wtcomb.npair/wtcomb.nlenses^2 $
                + wtrcomb.npair/wtrcomb.nlenses^2 )/wtcomb.area*1000.^2

  aploterror,!gratio,wtcomb.meanr[w],diffdensecomb[w],errorcomb[w],psym=4,$
    ytitle='Overdensity (# Mpc!U-2!N)'

  cumul = fltarr(nbin)
  FOR i=0L, nbin-1 DO BEGIN 
      IF i EQ 0 THEN BEGIN
          cumul[i] = (wtcomb.npair[i]/wtcomb.nlenses - wtrcomb.npair[i]/wtrcomb.nlenses)
      ENDIF ELSE BEGIN 
          cumul[i] = cumul[i-1] + (wtcomb.npair[i]/wtcomb.nlenses - wtrcomb.npair[i]/wtrcomb.nlenses)
      ENDELSE 
  ENDFOR 

  oplot,wtcomb.meanr[w], cumul
  legend,['Overdensity','Cumulative'],psym=[4,0],/center,/top,charsize=1.0,$
    thick=[!p.thick,!p.thick]


  !p.multi=[0,0,2]

  if dtype eq 'X' then key=get_kbrd(1)


  ;; fit gaussian

  tmptmp=gaussfit( wt.meanr[w[0:nbin-1]], diffdense[w[0:nbin-1]],a,nterms=4)
  tmptmp=gaussfit( wt2.meanr[w[0:nbin-1]], diffdense2[w[0:nbin-1]],a2,nterms=4)
  tmptmp=gaussfit( wtcomb.meanr[w[0:nbin-1]], diffdensecomb[w[0:nbin-1]],acomb,nterms=4)

  !p.multi=0

  nxx = 1000
  xx = arrscl( findgen(nxx), min(wt.meanr[w]), max(wt.meanr[w]) )

  z = (xx - a[1])/a[2]
  yy=a[0]*exp((-1.0*z^2)/2.0) + a[3]
  z = (xx - a2[1])/a2[2]
  yy2=a2[0]*exp((-1.0*z^2)/2.0) + a2[3]
  z = (xx - acomb[1])/acomb[2]
  yycomb=acomb[0]*exp((-1.0*z^2)/2.0) + acomb[3]

  aploterror,!gratio,wt.meanr[w],diffdense[w],error[w],psym=4,$
    title=stripestr[0],xtitle='radius (Mpc)',ytitle='Overdensity (#/Mpc!U2!N)',$
    yrange=[1.1*min(diffdense[w]),0.0]
  oplot,xx,yy
  if dtype eq 'X' then key=get_kbrd(1)
  aploterror,!gratio,wt2.meanr[w2],diffdense2[w2],error2[w2],psym=1,$
    title=stripestr[1],xtitle='radius (Mpc)',ytitle='Overdensity (#/Mpc!U2!N)',$
    yrange=[1.1*min(diffdense2[w]),0.0]
  oplot,xx,yy2
  if dtype eq 'X' then key=get_kbrd(1)
  aploterror,!gratio,wtcomb.meanr[w],diffdensecomb[w],errorcomb[w],psym=4,$
    xtitle='Projected Radius (h!U'+!tsym.minus+'1!N kpc)',$
    ytitle='Overdensity (# Mpc!U'+!tsym.minus+'2!N)',$
    yrange=[1.1*min(diffdensecomb[w]),0.0]
  oplot,xx,yycomb
  
  print,'stripe '+stripestr[0]
  print,a
  print,'stripe '+stripestr[1]
  print,a2
  print,'all'
  print,acomb

  meanr = wtcomb.meanr[w]
  diffdensecomb=diffdensecomb[w]
  errorcomb=errorcomb[w]
  cumul=cumul[w]


END 
