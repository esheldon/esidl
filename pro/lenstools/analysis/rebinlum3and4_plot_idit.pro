PRO rebinlum3and4_plot_idit

  lumclr=2

  combdir = $
    '/net/cheops2/home/esheldon/lensout/combstripe/comb/sublum/'+$
    !colors[lumclr]+'/'
  lumstr = 'fourbin'

  f1 = combdir + 'lum1'+lumstr+'_zgal_gal_stripe10_11_12_35_36_37_ri_jack_comb_N4.fit'
  f2 = combdir + 'lum2'+lumstr+'_zgal_gal_stripe10_11_12_35_36_37_ri_jack_comb_N4.fit'
  f3 = combdir + 'lum3'+lumstr+'_zgal_gal_stripe10_11_12_35_36_37_ri_jack_comb_N4.fit'
  f4 = combdir + 'lum4'+lumstr+'_zgal_gal_stripe10_11_12_35_36_37_ri_jack_comb_N4.fit'

  t1=mrdfits(f1,1) 
  t2=mrdfits(f2,1)
  t3=mrdfits(f3,1)
  t4=mrdfits(f4,1)

  combine_zshear, [f3, f4], tcomb

  zshear_rebin, tcomb, meanr_rebin, sigma_rebin, sigmaerr_rebin
  tcomb.sigma_rebin = sigma_rebin
  tcomb.sigmaerr_rebin = sigmaerr_rebin

  denscont2wgm, tcomb, wgm, wgmerr, r0, r0errl, r0errh, r0low, r0high
  denscont2wgm, tcomb, wgm_rebin, wgmerr_rebin, /rebin

  wgmfac = wgm[0]/tcomb.sigma[0]
  wgm1fac = t1.wgm[0]/t1.sigma[0]
  wgm2fac = t2.wgm[0]/t2.sigma[0]

  fitpower, tcomb.meanr/1000., wgm, wgmerr, [150., -0.8], tyfit, Aout, Asig

  aploterror, 1,tcomb.meanr/1000., wgm, wgmerr, /xlog,/ylog,$
              xrange=[0.01,20.0],/xsty,yrange=[1,20000],/yst,psym=8
  oplot,tcomb.meanr/1000.,tyfit
  
  IF !d.name EQ 'X' THEN key=get_kbrd(1)

  rangefac = 10.0; ELSE rangefac=5.0
  Asig[1] = Asig[1]*rangefac
  Asig[0] = Asig[0]*rangefac

  spowrange = Aout[1] + [-Asig[1], +Asig[1]]
  normrange = Aout[0] + [Asig[0], -Asig[0]]
  nnorm = 200
  npow = 200

  covariance = t1.covariance
  fac = (tcomb.tsigmaerr[17]/t1.tsigmaerr[17])^2
  covariance = covariance*fac*wgmfac^2

  pow_chisq_conf_gen, tcomb.meanr/1000., wgm, covariance, $
                      spowrange, normrange, npow, nnorm, $
                      chisq_surf, $
                      bestpow, bestnorm, $
                      powlow, powhigh, $
                      normlow, normhigh, $
                      yallow_low = yallow_low, $
                      yallow_high= yallow_high, yfit=yfit, /nodisplay

  aploterror, 1, tcomb.meanr_rebin/1000., wgm_rebin, wgmerr_rebin,$
              /xlog, /ylog, psym=8, $
              yrange = [1,20000], xrange=[10, 20000]/1000., /xst,/yst, $
              xtitle=!mpcxtitle,ytitle=!wgmytitle

  oploterror, tcomb.meanr_rebin/1000., wgm_rebin, wgmerr_rebin, psym=8, $
              color=!red,errc=!red
  oploterror, tcomb.meanr_rebin[6:8]/1000., wgm_rebin[6:8], wgmerr_rebin[6:8], psym=8, $
              color=!red,errc=!red
  oplot, tcomb.meanr/1000., yfit, color=!red

  oploterror, t1.meanr_rebin/1000., t1.wgm_rebin, t1.wgmerr_rebin, psym=8,$
              color=!blue,errc=!blue
  oploterror, t1.meanr_rebin[6:8]/1000., t1.wgm_rebin[6:8], t1.wgmerr_rebin[6:8], psym=8,$
              color=!blue,errc=!blue
  oplot, t1.meanr/1000., t1.yfit*wgm1fac, color=!blue

  oploterror, t2.meanr_rebin/1000., t2.wgm_rebin, t2.wgmerr_rebin, psym=8,$
              color=!green,errc=!green
  oploterror, t2.meanr_rebin[6:8]/1000., t2.wgm_rebin[6:8], t2.wgmerr_rebin[6:8], psym=8,$
              color=!green,errc=!green
  oplot, t2.meanr/1000., t2.yfit*wgm2fac, color=!green

  IF !d.name EQ 'X' THEN key=get_kbrd(1)

  aplot, 1, [0],/nodata, /xlog, /ylog, $
         yrange = [1,20000], xrange=[10, 20000]/1000., /xst,/yst, $
         xtitle=!mpcxtitle,ytitle=!wgmytitle

  x = [t1.meanr,reverse(t1.meanr)]/1000.
  y = [t1.yallow_low, reverse(t1.yallow_high)]*wgm1fac
  polyfill, x, y, color=!blue,/data

  x = [t2.meanr,reverse(t1.meanr)]/1000.
  y = [t2.yallow_low, reverse(t2.yallow_high)]*wgm2fac
  polyfill, x, y, color=!green,/data

  x = [tcomb.meanr,reverse(tcomb.meanr)]/1000.
  y = [yallow_low, reverse(yallow_high)]
  polyfill, x, y, color=!red,/data

  IF !d.name EQ 'X' THEN key=get_kbrd(1)

  ;; idit's stuff
  iabsmaglow = [-20., -21.5, -23.]
  iabsmaghigh = [-18.5, -20.0, -21.5]
  iabsmag=(iabsmaghigh+iabsmaglow)/2.
      
  ir0 = [4.72,6.28,7.42]
  ir0err = [0.44,0.77,0.33]

  r0 = [t1.r0, t2.r0, r0]
  r0low = [t1.r0low, t2.r0low, r0low]
  r0high = [t1.r0high, t2.r0high, r0high]

  wsum3 = total(t3.wsum)
  wsum4 = total(t4.wsum)
  wsum = wsum3 + wsum4
  meanlum = ( t3.tmeanlum*wsum3 + t4.tmeanlum*wsum4)/wsum
  meanlum = [t1.tmeanlum, t2.tmeanlum, meanlum]

  sun=[6.38,5.06,4.64,4.53,4.52]
  absmag = sun[lumclr] - 2.5*alog10(meanlum*1.e10)

  xrange = [-18, -24]
  xtitle = 'M!Dr!N - 5 log!D10!Nh'
  ytitle = 'r!D0!N [h!U'+!csym.minus+'1!N Mpc]'
  aplotderror, 1, iabsmag,ir0,iabsmaglow,iabsmaghigh,ir0-ir0err,ir0+ir0err,psym=8,$
             xrange=xrange, /xsty, yrange=[2,10],xtitle=xtitle,ytitle=ytitle

  IF !d.name EQ 'X' THEN cc = !green ELSE cc = !blue
  absmaglow = fltarr(3)
  absmaghigh = fltarr(3)

  diff = (absmag[1]-absmag[0])/2.
  absmaglow[0] = absmag[0]-diff
  absmaghigh[0] = absmag[0] + diff
  absmaglow[1] = absmag[1]-diff

  diff = (absmag[2]-absmag[1])/2.
  absmaghigh[1] = absmag[1] + diff
  absmaglow[2] = absmag[2] - diff
  absmaghigh[2] = absmag[2] + diff


  forprint,absmag,absmaglow,absmaghigh

  oplotderror, absmag, r0, absmaglow, absmaghigh, r0low,r0high, psym=8,$
               color=cc,errc=cc

END 
