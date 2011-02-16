PRO pow_scale_fit_one_bandpass, datax, data, dataerr, power, norm, $
                            _extra=extra

  pow_chisq_conf, datax, data, dataerr, power, norm, $
    chisq_surf, bestpow, bestnorm, powlow, powhigh, normlow, normhigh, $
    xtitle='!7a!X', ytitle='A', _extra=extra
  clevel = 0
  
  IF n_elements(bestpow) EQ 0 THEN return

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
    xtitle='!7a!X', ytitle='A', _extra=extra

  IF n_elements(bestpow) EQ 0 THEN return

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


;  rotate_plot,chisq_surf

  return
END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Main
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO fit_mass_vs_tag_nobin, glensum, rlensum, ilensum, tag, nrad, normrange, powrange, nnorm=nnorm,npower=npower,fiteach=fiteach

  IF n_params() LT 7 THEN BEGIN 
      print,'-Syntax: fit_mass_vs_tag_nobin, glensum, rlensum, ilensum, tag, nrad, normrange, powrange [, nnorm=nnorm,npower=npower,fiteach=fiteach]'
      return
  ENDIF 

  IF n_elements(nnorm) EQ 0 THEN nnorm=50L
  IF n_elements(npower) EQ 0 THEN npower=50L

  ;; need tag, maxtag, or xrange for model, normrange, powrange,  nx, npower, nnorm

  ;; this is to test fitting to scatter mass plot using
  ;; my chisq fitter (extremely slow for this large data set)

  ;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;

  IF NOT tag_exist(glensum, tag, index=wt) THEN BEGIN 
      print,'Tag '+tag+' not found'
      return
  ENDIF 

  colors = ['u','g','r','i','z']

  IF (!d.flags AND 1) EQ 0 THEN doX=1 ELSE doX=0
  IF doX THEN BEGIN
      !p.charsize = 1.0
      charsize = 1.2
      !p.thick = 1
      !x.thick = 1
      !y.thick = 1
      !p.charthick=1
  ENDIF ELSE BEGIN
      charsize=1.0
      !p.charsize = 1.0
      !p.thick = 5
      !x.thick = 5
      !y.thick = 5
      !p.charthick = 4
  ENDELSE 

  simpctable
  !p.background=!white
  !p.color = !black

  time=systime(1)

  ng = n_elements(glensum)
  nr = n_elements(rlensum)
  ni = n_elements(ilensum)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Match up the g,r,i for fitting combined
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  gindex = lindgen(ng)  
  rindex = lindgen(nr)  
  iindex = lindgen(ni)
  wg=where(glensum.sismasserr[nrad] gt 0.,ngrn)
  gindex = gindex[wg]
  wr=where(rlensum.sismasserr[nrad] gt 0.,nred)
  rindex = rindex[wr]
  wi=where(ilensum.sismasserr[nrad] gt 0.,nii)
  iindex = iindex[wi]

  print,'Starting with ',ngrn,nred,nii

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

  print
  print,ntostr(n_elements(iindex)),' matches found'
  print

  IF (n_elements(rindex) NE n_elements(iindex)) OR (n_elements(rindex) NE n_elements(gindex)) THEN BEGIN
      print,'Not same'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; get tag/mass values for each lens
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  g_tag = glensum.(wt)
  r_tag = rlensum.(wt)
  i_tag = ilensum.(wt)

  gmass = glensum.sismass[nrad]/1.e12
  rmass = rlensum.sismass[nrad]/1.e12
  imass = ilensum.sismass[nrad]/1.e12
  gmasserr = glensum.sismasserr[nrad]/1.e12
  rmasserr = rlensum.sismasserr[nrad]/1.e12
  imasserr = ilensum.sismasserr[nrad]/1.e12

  maxtag = max(g_tag) & mintag = min(g_tag)
  IF mintag LT 0. THEN xmin = 1.1*mintag ELSE xmin = 0.9*mintag
  IF maxtag LT 0. THEN xmax = 0.9*maxtag ELSE xmax = 1.1*maxtag

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; create the model
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  norm = arrscl( findgen(nnorm), normrange[0], normrange[1] )
  power = arrscl( findgen(npower), powrange[0], powrange[1] )
 
  IF keyword_set(fiteach) THEN BEGIN 
      print,'Fitting to g-band mass'
      pow_scale_fit_one_bandpass, g_tag[wg], gmass[wg], gmasserr[wg], power, norm, $
        title='g-band mass'
      IF doX THEN print,'Hit a key' & key=get_kbrd(1)
      print,'Fitting to r-band mass'
      pow_scale_fit_one_bandpass, r_tag[wr], rmass[wr], rmasserr[wr], power, norm, $
        title='r-band mass'
      IF doX THEN print,'Hit a key' & key=get_kbrd(1)
      print,'Fitting to i-band mass'
      pow_scale_fit_one_bandpass, i_tag[wi], imass[wi], imasserr[wi], power, norm, $
        title='i-band mass'
      IF doX THEN print,'Hit a key' & key=get_kbrd(1)
  ENDIF 
  
  time=systime(1)
  print,'Fitting to combined mass'
  pow_scale_fit_three_bandpass, g_tag[gindex], $
    gmass[gindex], gmasserr[gindex], $
    rmass[rindex], rmasserr[rindex], $
    imass[iindex], imasserr[iindex], $
    power, norm, $
    title='combined mass'
  ptime,systime(1)-time
      
  return

END 
