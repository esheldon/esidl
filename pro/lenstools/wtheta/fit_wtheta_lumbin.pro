PRO fit_wtheta_lumbin,clr,lumbin

  IF n_params() LT 2 THEN BEGIN
      print,'-Syntax: fit_wtheta_lumbin,clr,lumbin'
      return
  ENDIF 

  setupplot, dtype
  setup_mystuff

  lumbinstr='lum'+ntostr(long(lumbin))
  type=lumbinstr+'_wtheta'
  stripestr = ['stripe10','stripe82','stripe36','stripe37','stripe42','stripe43','all']
  nstripe = n_elements(stripestr) - 1

  indir='/sdss5/data0/lensout/'+stripestr+'/sublum/'+!COLORS[clr]+'/'
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

  cumul = fltarr(nbin)
  FOR i=0L, nbin-1 DO BEGIN 
      cumul[i] = wtcomb.tnpair[i]/wtcomb.nlenses-wtrcomb.tnpair[i]/wtrcomb.nlenses
  ENDFOR 
  aploterror,!gratio,wtcomb.meanr[w],diffdensecomb[w],errorcomb[w],psym=4,$
    ytitle='Overdensity (# Mpc!U-2!N)',xtitle='Projected Radius (kpc)',$
    yrange=[0,1.1*max(diffdensecomb[w])]
  oplot,wtcomb.rmax_act[w],cumul
  legend,['Overdensity','Cumulative'],psym=[4,0],/center,/top,charsize=1.0,$
    thick=[!p.thick,!p.thick],box=0


  !p.multi=[0,0,2]


  IF dtype EQ 'X' THEN key=get_kbrd(1)
  aguess = [1.0, -0.7]


  minpow = -1.5 & maxpow = -0.1
  minnorm = 0.1 & maxnorm = 1.5
  npow = 100L
  nnorm = 100L
  powvals = arrscl( findgen(npow), minpow, maxpow )
  normvals = arrscl( findgen(nnorm), minnorm, maxnorm )
  print

  pow_chisq_conf, wt.meanr[w[1:nbin-1]]/1000., $
    diffdense[w[1:nbin-1]], $
    error[w[1:nbin-1]], $
    powvals, normvals, $
    chisq_surf, $
    pmin, nmin, $
    powlow, powhigh, normlow, normhigh, $
    powallow=powallow,normallow=normallow
  print,'norm = ',nmin
  print,'pow = ',pmin

  pow_chisq_conf, wt2.meanr[w2[1:nbin2-1]]/1000., $
    diffdense2[w2[1:nbin2-1]], $
    error2[w2[1:nbin2-1]], $
    powvals, normvals, $
    chisq_surf2, $
    pmin2, nmin2, $
    powlow2, powhigh2, normlow2, normhigh2, $
    powallow=powallow2,normallow=normallow2
  print,'norm2 = ',nmin2
  print,'pow2 = ',pmin2

  if dtype eq 'X' then key=get_kbrd(1)
  pow_chisq_conf, wtcomb.meanr[w[1:nbin-1]]/1000., $
    diffdensecomb[w[1:nbin-1]], $
    errorcomb[w[1:nbin-1]], $
    powvals, normvals, $
    chisq_surfcomb, $
    pmincomb, nmincomb, $
    powlowcomb, powhighcomb, normlowcomb, normhighcomb, $
    powallow=powallowcomb,normallow=normallowcomb
  print,'norm comb = ',nmincomb
  print,'pow comb = ',pmincomb

  !p.multi=0

  nxx = 1000
  xx = arrscl( findgen(nxx), min(wt.meanr[w])/1000., max(wt.meanr[w])/1000. )
  yy = nmin*(xx)^pmin
  yy2 = nmin2*(xx)^pmin2
  yycomb = nmincomb*(xx)^pmincomb

  if dtype eq 'X' then key=get_kbrd(1)
  aploterror,!gratio,wt.meanr[w]/1000.,diffdense[w],error[w],psym=4,$
    title=stripestr[0],xtitle='radius (Mpc)',ytitle='Overdensity (#/Mpc!U2!N)',$
    yrange=[0,1.1*max(diffdense[w])]
  oplot,xx,yy
  if dtype eq 'X' then key=get_kbrd(1)
  aploterror,!gratio,wt2.meanr[w2]/1000.,diffdense2[w2],error2[w2],psym=1,$
    title=stripestr[1],xtitle='radius (Mpc)',ytitle='Overdensity (#/Mpc!U2!N)',$
    yrange=[0,1.1*max(diffdense2[w])]
  oplot,xx,yy2
  if dtype eq 'X' then key=get_kbrd(1)
  aploterror,!gratio,wtcomb.meanr[w]/1000.,diffdensecomb[w],errorcomb[w],psym=4,$
    title=stripestr[nstripe],xtitle='radius (Mpc)',ytitle='Overdensity (#/Mpc!U2!N)',$
    yrange=[0,1.1*max(diffdensecomb[w])]
  oplot,xx,yycomb

  if dtype eq 'X' then key=get_kbrd(1)

  modely_low = fltarr(nxx)
  modely_high = fltarr(nxx)
  modely_low2 = fltarr(nxx)
  modely_high2 = fltarr(nxx)
  modely_lowcomb = fltarr(nxx)
  modely_highcomb = fltarr(nxx)
  FOR ix=0L, nxx-1 DO BEGIN 

      mody = normallow*(xx[ix])^powallow
      modely_low[ix] = min(mody)
      modely_high[ix] = max(mody)

      mody = normallow2*(xx[ix])^powallow2
      modely_low2[ix] = min(mody)
      modely_high2[ix] = max(mody)

      mody = normallowcomb*(xx[ix])^powallowcomb
      modely_lowcomb[ix] = min(mody)
      modely_highcomb[ix] = max(mody)

  ENDFOR 

  lines = [0,0,0]

  xrange=[0.0,1.0]
  yrange=[0.0,8.0]
  aplot,!gratio,[0],/nodata,yrange=yrange,xrange=xrange

  FOR ix=0L, nxx-1 DO BEGIN 
      oplot,[ xx[ix], xx[ix] ], $
            [modely_low[ix], modely_high[ix]],line=lines[0]
      oplot,[ xx[ix], xx[ix] ], $
            [modely_low2[ix], modely_high2[ix]],line=lines[1],color=!red
      oplot,[ xx[ix], xx[ix] ], $
            [modely_lowcomb[ix], modely_highcomb[ix]],line=lines[2],color=!blue
  ENDFOR 
  ;oplot,xx, modely_high,line=0
  ;oplot,xx, modely_low,line=0
  ;oplot,xx, modely_high2, color=!red,line=3
  ;oplot,xx, modely_low2, color=!red,line=3
  ;oplot,xx, modely_highcomb,color=!yellow,line=2
  ;oplot,xx, modely_lowcomb,color=!yellow,line=2

  legend, [stripestr[0],stripestr[1],stripestr[nstripe]], line=lines,color=[!black, !red, !blue],/right,thick=[!p.thick,!p.thick,!p.thick]

  if dtype eq 'X' then key=get_kbrd(1)

  aploterror,!gratio,wtcomb.meanr[w],diffdensecomb[w],errorcomb[w],psym=4,$
    ytitle='Overdensity (# Mpc!U'+!tsym.minus+'2!N)',$
    xtitle='Projected Radius (h!U'+!tsym.minus+'1!N kpc)',$
    yrange=[0,1.1*max(diffdensecomb[w])]
  oplot,xx*1000.,yycomb

  oplot,wtcomb.rmax_act[w],cumul,line=2
  legend,['Overdensity','Cumulative'],line=[-1,2],/center,/top,charsize=1.0,$
    thick=[!p.thick,!p.thick],box=0
  legend,['Overdensity','Cumulative'],psym=[4,3],/center,/top,charsize=1.0,$
    thick=[!p.thick,!p.thick],box=0

END 
