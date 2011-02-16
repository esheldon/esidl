PRO scale_fit_one_bandpass, datax, data, dataerr, modelx, model, power, norm, $
                            _extra=extra

  chisq_conf, datax, data, dataerr, modelx, model, power, norm, $
    chisq_surf, bestpow, bestnorm, powlow, powhigh, normlow, normhigh, $
    xtitle='!7a!X', ytitle='A', _extra=extra
  clevel = 0
  
  highdiff = normhigh[clevel]-bestnorm & lowdiff = bestnorm - normlow[clevel]
  mess = 'A = '+ntostr(bestnorm,6)+$
    '!S!U+'+ntostr(highdiff,5)+'!R!D-'+ntostr(lowdiff,5)
  mess = [mess, '']
  highdiff = powhigh[clevel]-bestpow & lowdiff = bestpow - powlow[clevel]
  mess = [mess, $
          '!7a!X = '+ntostr(bestpow, 6) + $
          '!S!U+'+ntostr(highdiff,5)+'!R!D-'+ntostr(lowdiff,5)]
  legend, mess, /right, bcolor=bcolor, charsize=charsize
  
  print,'Norm: '+ntostr(bestnorm),' + ',$
    ntostr(normhigh[clevel]-bestnorm),' - ',ntostr(bestnorm-normlow[clevel])
  print,'Power: '+ntostr(bestpow),' + ',$
    ntostr(powhigh[clevel]-bestpow),' - ',ntostr(bestpow-powlow[clevel])
  return
END 

PRO scale_fit_three_bandpass, datax, $
                              gdata, gdataerr, rdata, rdataerr, idata, idataerr, $
                              modelx, model, power, norm, $
                              _extra=extra

  gg = 3.50959599
  rr = 2.235299689
  ii = 2.549197961
  gr = 1.586466533 
  gi = 1.6983151 
  ri = 1.512359722

  chisq_conf_3band, $
    datax, gdata, gdataerr, $
    datax, rdata, rdataerr, $
    datax, idata, idataerr, $
    gg, rr, ii, gr, gi, ri, $
    modelx, model, power, norm, $
    chisq_surf, bestpow, bestnorm, powlow, powhigh, normlow, normhigh, $
    xtitle='!7a!X', ytitle='A', _extra=extra

  clevel = 0
  highdiff = normhigh[clevel]-bestnorm & lowdiff = bestnorm - normlow[clevel]
  mess = 'A = '+ntostr(bestnorm,6)+$
    '!S!U+'+ntostr(highdiff,5)+'!R!D-'+ntostr(lowdiff,5)
  mess = [mess, '']
  highdiff = powhigh[clevel]-bestpow & lowdiff = bestpow - powlow[clevel]
  mess = [mess, $
          '!7a!X = '+ntostr(bestpow, 6) + $
          '!S!U+'+ntostr(highdiff,5)+'!R!D-'+ntostr(lowdiff,5)]
  legend, mess, /right, bcolor=bcolor, charsize=charsize
  
  print,'Norm: '+ntostr(bestnorm),' + ',$
    ntostr(normhigh[clevel]-bestnorm),' - ',ntostr(bestnorm-normlow[clevel])
  print,'Power: '+ntostr(bestpow),' + ',$
    ntostr(powhigh[clevel]-bestpow),' - ',ntostr(bestpow-powlow[clevel])
  return
END 

PRO pow_scale_fit_one_bandpass, datax, data, dataerr, power, norm, $
                                powallow=powallow,normallow=normallow, $
                                _extra=extra

  pow_chisq_conf, datax, data, dataerr, power, norm, $
    chisq_surf, bestpow, bestnorm, powlow, powhigh, normlow, normhigh, $
    powallow=powallow,normallow=normallow, $
    xtitle='!7a!X', ytitle='A', _extra=extra
  clevel = 0
  
  highdiff = normhigh[clevel]-bestnorm & lowdiff = bestnorm - normlow[clevel]
  mess = 'A = '+ntostr(bestnorm,6)+$
    '!S!U+'+ntostr(highdiff,5)+'!R!D-'+ntostr(lowdiff,5)
  mess = [mess, '']
  highdiff = powhigh[clevel]-bestpow & lowdiff = bestpow - powlow[clevel]
  mess = [mess, $
          '!7a!X = '+ntostr(bestpow, 6) + $
          '!S!U+'+ntostr(highdiff,5)+'!R!D-'+ntostr(lowdiff,5)]
  legend, mess, /right, bcolor=bcolor, charsize=charsize
  
  print,'Norm: '+ntostr(bestnorm),' + ',$
    ntostr(normhigh[clevel]-bestnorm),' - ',ntostr(bestnorm-normlow[clevel])
  print,'Power: '+ntostr(bestpow),' + ',$
    ntostr(powhigh[clevel]-bestpow),' - ',ntostr(bestpow-powlow[clevel])
  return
END 

PRO pow_scale_fit_three_bandpass, datax, $
                                  gdata, gdataerr, rdata, rdataerr, idata, idataerr, $
                                  power, norm, $
                                  powallow=powallow,normallow=normallow, $
                                  _extra=extra

  gg = 3.50959599
  rr = 2.235299689
  ii = 2.549197961
  gr = 1.586466533 
  gi = 1.6983151 
  ri = 1.512359722

  pow_chisq_conf_3band, $
    datax, gdata, gdataerr, $
    datax, rdata, rdataerr, $
    datax, idata, idataerr, $
    gg, rr, ii, gr, gi, ri, $
    power, norm, $
    chisq_surf, bestpow, bestnorm, powlow, powhigh, normlow, normhigh, $
    powallow=powallow,normallow=normallow, $
    xtitle='!7a!X', ytitle='A', _extra=extra

  clevel = 0
  highdiff = normhigh[clevel]-bestnorm & lowdiff = bestnorm - normlow[clevel]
  mess = 'A = '+ntostr(bestnorm,6)+$
    '!S!U+'+ntostr(highdiff,5)+'!R!D-'+ntostr(lowdiff,5)
  mess = [mess, '']
  highdiff = powhigh[clevel]-bestpow & lowdiff = bestpow - powlow[clevel]
  mess = [mess, $
          '!7a!X = '+ntostr(bestpow, 6) + $
          '!S!U+'+ntostr(highdiff,5)+'!R!D-'+ntostr(lowdiff,5)]
  legend, mess, /right, bcolor=bcolor, charsize=charsize
  
  print,'Norm: '+ntostr(bestnorm),' + ',$
    ntostr(normhigh[clevel]-bestnorm),' - ',ntostr(bestnorm-normlow[clevel])
  print,'Power: '+ntostr(bestpow),' + ',$
    ntostr(powhigh[clevel]-bestpow),' - ',ntostr(bestpow-powlow[clevel])
  return
END 

PRO fit_mass_vs_lum_nobin, glensum, rlensum, ilensum, nrad, create=create, each=each

  IF (n_params() LT 4) AND (NOT keyword_set(create)) THEN BEGIN 
      print,'-Syntax: fit_mass_vs_lum_nobin, glensum, rlensum, ilensum, nrad, create=create, each=each'
      return
  ENDIF 

  ;; this is to test fitting to scatter mass plot using
  ;; my chisq fitter 

  ;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;

  indir = '/sdss4/data1/esheldon/GAL_GAL/spectra/'
  IF keyword_set(create) THEN BEGIN 
      
      concat_structs,mrdfits(indir+'main_zgal_gal_stripe10_g_lensum_massadd_N1.fit',1),$
                     mrdfits(indir+'main_zgal_gal_stripe82_g_lensum_massadd_N1.fit',1),$
                     glensum
      concat_structs,mrdfits(indir+'main_zgal_gal_stripe10_r_lensum_massadd_N1.fit',1),$
                     mrdfits(indir+'main_zgal_gal_stripe82_r_lensum_massadd_N1.fit',1),$
                     rlensum
      concat_structs,mrdfits(indir+'main_zgal_gal_stripe10_i_lensum_massadd_N1.fit',1),$
                     mrdfits(indir+'main_zgal_gal_stripe82_i_lensum_massadd_N1.fit',1),$
                     ilensum
      return
  ENDIF 
  colors = ['u','g','r','i','z']

  IF (!d.flags AND 1) EQ 0 THEN doX=1 ELSE doX=0
  IF doX THEN BEGIN
      setupplot,'X'
  ENDIF ELSE BEGIN
      setupplot,'PS'
  ENDELSE 

  simpctable
  !p.background=!white
  !p.color = !black

  time=systime(1)
  wlum = [0,1,2,3,4]
;  wlum=3
  nlum = n_elements(wlum)
  
  ng = n_elements(glensum)
  nr = n_elements(rlensum)
  ni = n_elements(ilensum)

  maxlum = [10.0e10, 10.e10, 15.0e10, 30.0e10, 45.0e10]

  FOR ilum=0L, nlum-1 DO BEGIN 

      print
      print,'ng = ',ntostr(ng),' nr = ',ntostr(nr),' ni = ',ntostr(ni)
      print

      gindex = lindgen(ng)  
      rindex = lindgen(nr)  
      iindex = lindgen(ni)
      w=where( (glensum.sismasserr[nrad] GT 0.) AND $
               (glensum.lum[wlum[ilum]]  GT 0.) AND $
               (glensum.lum[wlum[ilum]]  LT maxlum[wlum[ilum]] ) )
      gindex = gindex[w]
      w=where( (rlensum.sismasserr[nrad] GT 0.) AND $
               (rlensum.lum[wlum[ilum]]  GT 0.) AND $
               (rlensum.lum[wlum[ilum]]  LT maxlum[wlum[ilum]] ) )
      rindex = rindex[w]
      w=where( (ilensum.sismasserr[nrad] GT 0.) AND $
               (ilensum.lum[wlum[ilum]]  GT 0.) AND $
               (ilensum.lum[wlum[ilum]]  LT maxlum[wlum[ilum]] ) )
      iindex = iindex[w]

      plothist, rlensum[rindex].lum[2]/1.e10,xhist1,yhist1,bin=0.5,line=2

      ;; match up
      photo_match, glensum[gindex].run,    glensum[gindex].rerun, $
                   glensum[gindex].camcol, glensum[gindex].field, $
                   glensum[gindex].id, $
                   rlensum[rindex].run,    rlensum[rindex].rerun, $
                   rlensum[rindex].camcol, rlensum[rindex].field, $
                   rlensum[rindex].id, $
                   gmatch, rmatch
      gindex = gindex[gmatch]
      rindex = rindex[rmatch]

      photo_match, glensum[gindex].run,    glensum[gindex].rerun, $
                   glensum[gindex].camcol, glensum[gindex].field, $
                   glensum[gindex].id, $
                   ilensum[iindex].run,    ilensum[iindex].rerun, $
                   ilensum[iindex].camcol, ilensum[iindex].field, $
                   ilensum[iindex].id, $
                   gmatch, imatch
      gindex = gindex[gmatch]
      iindex = iindex[imatch]

      photo_match, glensum[gindex].run,    glensum[gindex].rerun, $
                   glensum[gindex].camcol, glensum[gindex].field, $
                   glensum[gindex].id, $
                   rlensum[rindex].run,    rlensum[rindex].rerun, $
                   rlensum[rindex].camcol, rlensum[rindex].field, $
                   rlensum[rindex].id, $
                   gmatch, rmatch
      rindex = rindex[rmatch]

      print,ntostr(n_elements(iindex)),' matches found'
      print

      plothist, rlensum[rindex].lum[2]/1.e10,xhist2,yhist2,bin=0.5,/overplot
      plot,xhist1,yhist1-yhist2,line=10

return
      IF (n_elements(rindex) NE n_elements(iindex)) OR (n_elements(rindex) NE n_elements(gindex)) THEN BEGIN
          print,'Not same'
          return
      ENDIF 

      g_lum = glensum[gindex].lum[wlum[ilum]]/1.e10
      r_lum = rlensum[rindex].lum[wlum[ilum]]/1.e10
      i_lum = ilensum[iindex].lum[wlum[ilum]]/1.e10
      gmass = glensum[gindex].sismass[nrad]/1.e12
      rmass = rlensum[rindex].sismass[nrad]/1.e12
      imass = ilensum[iindex].sismass[nrad]/1.e12
      gmasserr = glensum[gindex].sismasserr[nrad]/1.e12
      rmasserr = rlensum[rindex].sismasserr[nrad]/1.e12
      imasserr = ilensum[iindex].sismasserr[nrad]/1.e12

      nnorm  = 75L
      npower = 75L
      CASE wlum[ilum] OF 
          0:  powrange = [-0.5, 1.5]
          1:  powrange = [0.0, 2.5]
          2:  powrange = [0.0, 2.5]
          3:  powrange = [0.0, 2.5]
          4:  powrange = [0.0, 2.5]
      ENDCASE
      rmax_3rdbin = mean(rlensum.rmax_act[2])
      rmax_thisbin = mean(rlensum.rmax_act[nrad])
      norm_fac = rmax_thisbin/rmax_3rdbin
      CASE wlum[ilum] OF 
          0:  normrange = [0.0,6.0]*norm_fac
          1:  normrange = [0.0,6.0]*norm_fac
          2:  normrange = [0.,5.]*norm_fac
          3:  normrange = [0.,4.]*norm_fac
          4:  normrange = [0.,4.]*norm_fac
      ENDCASE
 
      norm = arrscl( findgen(nnorm), normrange[0], normrange[1] )
      power = arrscl( findgen(npower), powrange[0], powrange[1] )

      ;; new power-only fitter
      IF keyword_set(each) THEN BEGIN 
          time=systime(1)
          print,'Fitting to g-band mass    '+colors[wlum[ilum]]+'-band luminosity'
          pow_scale_fit_one_bandpass, g_lum, gmass, gmasserr, power, norm, $
            title='g-band mass   '+colors[wlum[ilum]]+'-band luminosity'
          ptime,systime(1)-time

          time=systime(1)
          print,'Fitting to r-band mass    '+colors[wlum[ilum]]+'-band luminosity'
          pow_scale_fit_one_bandpass, r_lum, rmass, rmasserr, power, norm, $
            title='r-band mass   '+colors[wlum[ilum]]+'-band luminosity'
          ptime,systime(1)-time

          time=systime(1)
          print,'Fitting to i-band mass    '+colors[wlum[ilum]]+'-band luminosity'
          pow_scale_fit_one_bandpass, i_lum, imass, imasserr, power, norm, $
            title='i-band mass   '+colors[wlum[ilum]]+'-band luminosity'
          ptime,systime(1)-time
      ENDIF 
      time=systime(1)
      print,'Fitting to combined mass   '+colors[wlum[ilum]]+'-band luminosity'
      pow_scale_fit_three_bandpass, g_lum, $
        gmass, gmasserr, rmass, rmasserr, imass, imasserr, $
        power, norm, $
        title='combined mass   '+colors[wlum[ilum]]+'-band luminosity', $
        powallow=powallow,normallow=normallow
      ptime,systime(1)-time

      ;; write out the allowed region
      outfile = indir+'massvs_'+colors[wlum[ilum]]+'lum_allowedpar.fit'
      tmpstruct = create_struct('powallow',powallow,$
                                'normallow',normallow)
      print,'Allowed parameter file: ',outfile
      mwrfits, tmpstruct, outfile, /create

  ENDFOR 
  return

END 
