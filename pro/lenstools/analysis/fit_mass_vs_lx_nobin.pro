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

PRO fit_mass_vs_lx_nobin, nrad, glensum, rlensum, ilensum, create=create, lxtag=lxtag, maxz=maxz

  IF (n_params() LT 4) AND (NOT keyword_set(create)) THEN BEGIN 
      print,'-Syntax: test, nrad, glensum, rlensum, ilensum, create=create'
      return
  ENDIF 

  ;; this is to test fitting to scatter mass plot using
  ;; my chisq fitter (extremely slow for this large data set)
  ;; this is a keeper

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

  IF n_elements(maxz) EQ 0 THEN maxz = 0.4

  IF n_elements(lxtag) EQ 0 THEN lxtag = 'l_x_bol'
  IF NOT tag_exist(glensum,lxtag,index=wlx) THEN BEGIN 
      print,'No such tag'
      return
  ENDIF 

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

  stripes = ['752_756','94_125']

  nstripe = n_elements(stripes)
  clr = [1,2,3]
  colors = ['u','g','r','i','z']
  nclr = n_elements(clr)

  dens = fltarr(nstripe, nclr)
  dens[0,0] = 1.2312680         ;g
  dens[0,1] = 1.7684189         ;r
  dens[0,2] = 1.5075363         ;i
  dens[1,0] = 0.80527492        ;g
  dens[1,1] = 1.3028618         ;r
  dens[1,2] = 1.0964445         ;i
  stripes = ['752_756']
  nstripe=1
  IF n_elements(glensum) EQ 0 THEN BEGIN 
;      indir = '/sdss4/data1/esheldon/CLUSTER/RASS/BOOT/chris/newcat/'
      indir = '/sdss4/data1/esheldon/CLUSTER/RASS/BOOT/lambda_cosmo/rmax2000/'
      
      ginf = indir+'subz_RASS_'+stripes+'_g_lensum_N1.fit'
      rinf = indir+'subz_RASS_'+stripes+'_r_lensum_N1.fit'
      iinf = indir+'subz_RASS_'+stripes+'_i_lensum_N1.fit'
              
      goutf = indir+'subz_RASS_'+stripes+'_g_lensum_massadd_N1.fit'
      routf = indir+'subz_RASS_'+stripes+'_r_lensum_massadd_N1.fit'
      ioutf = indir+'subz_RASS_'+stripes+'_i_lensum_massadd_N1.fit'

      IF keyword_set(create) THEN BEGIN 
          FOR ist=0L, nstripe-1 DO BEGIN 
              
              ginfile = ginf[ist]
              rinfile = rinf[ist]
              iinfile = iinf[ist]
              
              goutfile = goutf[ist]
              routfile = routf[ist]
              ioutfile = ioutf[ist]

              glsum = mrdfits(ginfile,1,ghdr)
              rlsum = mrdfits(rinfile,1,rhdr)
              ilsum = mrdfits(iinfile,1,ihdr)

              fxhclean, ghdr
              fxhclean, rhdr
              fxhclean, ihdr

              binsize = sxpar(ghdr, 'binwidth') & print,binsize
              rminkpc = sxpar(ghdr, 'rminkpc') & print,rminkpc
              rmaxkpc = sxpar(ghdr, 'rmaxkpc') > 1000. & print,rmaxkpc
              hval = sxpar(ghdr,'h')
              print,'Adding mass to glensum  ',stripes[ist]
              lensum_sis_fit, temporary(glsum), dens[ist,0], binsize, $
                rminkpc, rmaxkpc, hval, glensum
              print,'Adding mass to rlensum  ',stripes[ist]
              lensum_sis_fit, temporary(rlsum), dens[ist,1], binsize, $
                rminkpc, rmaxkpc, hval, rlensum
              print,'Adding mass to ilensum  ',stripes[ist]
              lensum_sis_fit, temporary(ilsum), dens[ist,2], binsize, $
                rminkpc, rmaxkpc, hval, ilensum
              
              mwrfits, temporary(glensum), goutfile, ghdr, /create
              mwrfits, temporary(rlensum), routfile, rhdr, /create
              mwrfits, temporary(ilensum), ioutfile, ihdr, /create
          ENDFOR 
          return
      ENDIF ELSE BEGIN 

          FOR ist=0L, nstripe-1 DO BEGIN 
                            
              goutfile = goutf[ist]
              routfile = routf[ist]
              ioutfile = ioutf[ist]

              glsum = mrdfits(goutfile,1,ghdr)
              rlsum = mrdfits(routfile,1,rhdr)
              ilsum = mrdfits(ioutfile,1,ihdr)


              IF ist EQ 0 THEN BEGIN 
                  glensum = temporary(glsum)
                  rlensum = temporary(rlsum)
                  ilensum = temporary(ilsum)
              ENDIF ELSE BEGIN 
                  concat_dstructs, temporary(glsum),temporary(glensum),tmp
                  glensum = temporary(tmp)
                  concat_dstructs, temporary(rlsum),temporary(rlensum),tmp
                  rlensum = temporary(tmp)
                  concat_dstructs, temporary(ilsum),temporary(ilensum),tmp
                  ilensum = temporary(tmp)
              ENDELSE 
          ENDFOR 
      ENDELSE 
  ENDIF 

  ng = n_elements(glensum)
  nr = n_elements(rlensum)
  ni = n_elements(ilensum)

  gindex = lindgen(ng)  
  rindex = lindgen(nr)  
  iindex = lindgen(ni)
  w=where(glensum.sismasserr[nrad] gt 0. AND glensum.z LT maxz,ngrn)
  gindex = gindex[w]
  w=where(rlensum.sismasserr[nrad] gt 0. AND rlensum.z LT maxz,nred)
  rindex = rindex[w]
  w=where(ilensum.sismasserr[nrad] gt 0. AND ilensum.z LT maxz,nii)
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

  g_lx = glensum[gindex].(wlx)
  r_lx = rlensum[rindex].(wlx)
  i_lx = ilensum[iindex].(wlx)
  gmass = glensum[gindex].sismass[nrad]/1.e14
  rmass = rlensum[rindex].sismass[nrad]/1.e14
  imass = ilensum[iindex].sismass[nrad]/1.e14
  gmasserr = glensum[gindex].sismasserr[nrad]/1.e14
  rmasserr = rlensum[rindex].sismasserr[nrad]/1.e14
  imasserr = ilensum[iindex].sismasserr[nrad]/1.e14
      
;  fitlin, r_lx, rmass, rmasserr

  nx = 500.
  nnorm  = 50L
  npower = 50L

  powrange = [-1., 2.0]
  normrange = [0., 12.]

  xmin = 0.001
  xmax = 1.1*max(g_lx)

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
  scale_fit_one_bandpass, g_lx, gmass, gmasserr, modelx, model, power, norm, $
    title='g-band mass'
  print,'Fitting to r-band mass'
  scale_fit_one_bandpass, r_lx, rmass, rmasserr, modelx, model, power, norm, $
    title='r-band mass'
  print,'Fitting to i-band mass'
  scale_fit_one_bandpass, i_lx, imass, imasserr, modelx, model, power, norm, $
    title='i-band mass'

  print,'Fitting to combined mass'
  scale_fit_three_bandpass, g_lx, $
    gmass, gmasserr, rmass, rmasserr, imass, imasserr, $
    modelx, model, power, norm, $
    title='combined mass'

      
  return

END 
