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

PRO fit_mass_vs_dn_nobin, nrad, glensum, rlensum, ilensum, create=create

  IF (n_params() LT 4) AND (NOT keyword_set(create)) THEN BEGIN 
      print,'-Syntax: test, nrad, glensum, rlensum, ilensum, create=create'
      return
  ENDIF 

  ;; this is to test fitting to scatter mass plot using
  ;; my chisq fitter (extremely slow for this large data set)

  ;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;

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
  wlum = [0,1,2,3,4]
;  wlum=3
  nlum = n_elements(wlum)

  gdens = 1.2312680             ;ngal/square arcminute
  rdens = 1.7684189       
  idens = 1.5075363
  
  IF n_elements(glensum) EQ 0 THEN BEGIN 
      indir = '/sdss4/data1/esheldon/GAL_GAL/spectra/'
      ginfile = indir+'main_zgal_gal_752_756_g_lensum_massadd_N1.fit'
      rinfile = indir+'main_zgal_gal_752_756_r_lensum_massadd_N1.fit'
      iinfile = indir+'main_zgal_gal_752_756_i_lensum_massadd_N1.fit'

      glensum = mrdfits(ginfile, 1, ghdr)
      rlensum = mrdfits(rinfile, 1, rhdr)
      ilensum = mrdfits(iinfile, 1, ihdr)
      fit_dn,glensum,dn,clr=2.,thresh=22.
      glensum.dn = dn
      fit_dn,rlensum,dn,clr=2.,thresh=22.
      rlensum.dn = dn
      fit_dn,ilensum,dn,clr=2.,thresh=22.
      ilensum.dn = dn
  ENDIF 

  ng = n_elements(glensum)
  nr = n_elements(rlensum)
  ni = n_elements(ilensum)

  nlum = 1                      ;only one dn
  maxdn = 10.
  FOR ilum=0L, nlum-1 DO BEGIN 

      gindex = lindgen(ng)  
      rindex = lindgen(nr)  
      iindex = lindgen(ni)
      w=where(glensum.sismasserr[nrad] gt 0. and glensum.dn GT 0. AND glensum.dn LT maxdn,ngrn)
      gindex = gindex[w]
      w=where(rlensum.sismasserr[nrad] gt 0. and rlensum.dn GT 0. AND rlensum.dn LT maxdn,nred)
      rindex = rindex[w]
      w=where(ilensum.sismasserr[nrad] gt 0. and ilensum.dn GT 0. AND ilensum.dn LT maxdn,nii)
      iindex = iindex[w]

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

;      IF (n_elements(rindex) NE n_elements(iindex)) OR (n_elements(rindex) NE n_elements(gindex)) THEN BEGIN
;          print,'Not same'
;          return
;      ENDIF 

      g_dn = glensum[gindex].dn
      r_dn = rlensum[rindex].dn
      i_dn = ilensum[iindex].dn
      gmass = glensum[gindex].sismass[nrad]/1.e12
      rmass = rlensum[rindex].sismass[nrad]/1.e12
      imass = ilensum[iindex].sismass[nrad]/1.e12
      gmasserr = glensum[gindex].sismasserr[nrad]/1.e12
      rmasserr = rlensum[rindex].sismasserr[nrad]/1.e12
      imasserr = ilensum[iindex].sismasserr[nrad]/1.e12
      
;      help,g_lum,r_lum
;      w=where(g_lum NE r_lum,nw)
;      print,nw
;      print,max(g_lum - r_lum),max(g_lum-i_lum)
;      return

      nx = 500.
      nnorm  = 30L
      npower = 30L

      powrange = [0.0, 2.0]
      normrange = [0., 6.]

      xmin = 0.
      xmax = maxdn+0.5          ;kpc
;      print,'Fit g'
;      fitlin, g_lum, gmass, gmasserr
;      print,'Fit r'
;      fitlin, r_lum, rmass, rmasserr
;      print,'Fit i'
;      fitlin, i_lum, imass, imasserr
     

      IF n_elements(oldrange) EQ 0 THEN BEGIN
          oldrange=normrange
          powmodel, normrange, powrange, xmin, xmax, modelx, model, norm, power,$
            nnorm=nnorm,npower=npower,nx=nx
      ENDIF  ELSE BEGIN 
          IF (normrange[0] NE oldrange[0]) OR (normrange[1] NE oldrange[1]) THEN BEGIN
              powmodel, normrange, powrange, xmin, xmax, modelx, model, norm, power,$
                nnorm=nnorm,npower=npower,nx=nx
          ENDIF 
      ENDELSE 
      oldrange = normrange

      print,'Fitting to g-band mass'
      scale_fit_one_bandpass, g_dn, gmass, gmasserr, modelx, model, power, norm, $
        title='g-band mass'
      print,'Fitting to r-band mass'
      scale_fit_one_bandpass, r_dn, rmass, rmasserr, modelx, model, power, norm, $
        title='r-band mass'
      print,'Fitting to i-band mass'
      scale_fit_one_bandpass, i_dn, imass, imasserr, modelx, model, power, norm, $
        title='i-band mass'

      print,'Fitting to combined mass'
      scale_fit_three_bandpass, g_dn, $
        gmass, gmasserr, rmass, rmasserr, imass, imasserr, $
        modelx, model, power, norm, $
        title='combined mass'

      
  ENDFOR 
  return

END 
