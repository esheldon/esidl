PRO make_lrg_photoz_plots, scat, lcat, nzstruct, w, nostar=nostar

  ;; testing lrg selection

  stripes = [9,10,11,12,13,14,15]
  clrs = [1,2,3]

  IF n_elements(scat) EQ 0 THEN BEGIN 
      get_scat,stripes,clrs,scat, /hirata
  ENDIF 

  IF n_elements(lcat) EQ 0 THEN BEGIN 
      rdlensout, lcat
  ENDIF 

  IF n_elements(nzstruct) EQ 0 THEN BEGIN 
      get_nz_photoz, stripes, clrs, nzstruct, /hirata
  ENDIF 
  wsource = nzstruct.useind
  make_lrg_cuts, scat, lrg_flags, nostar=nostar

  ;; Select lrg of various redshifts
  make_lrg_struct, lrgs
  lrgs.lrg = 'Y'
  lrg_select, blah, lrgs, w_lrg, $
    lrg_flags=lrg_flags, input_index=wsource
  n_lrg = n_elements(w_lrg)

  make_lrg_struct, lrgs
  lrgs.lrg_lo_med = 'Y'
  lrg_select, blah, lrgs, w_lrg_lo_med, $
    lrg_flags=lrg_flags, input_index=wsource
  n_lrg_lo_med = n_elements(w_lrg_lo_med)

  make_lrg_struct, lrgs
  lrgs.lrg_hi_med = 'Y'
  lrg_select, blah, lrgs, w_lrg_hi_med, $
    lrg_flags=lrg_flags, input_index=wsource  
  n_lrg_hi_med = n_elements(w_lrg_hi_med)

  make_lrg_struct, lrgs
  lrgs.lrg_lo_deep = 'Y'
  lrg_select, blah, lrgs, w_lrg_lo_deep, $
    lrg_flags=lrg_flags, input_index=wsource
  n_lrg_lo_deep = n_elements(w_lrg_lo_deep)

  make_lrg_struct, lrgs
  lrgs.lrg_hi_deep = 'Y'
  lrg_select, blah, lrgs, w_lrg_hi_deep, $
    lrg_flags=lrg_flags, input_index=wsource
  n_lrg_hi_deep = n_elements(w_lrg_hi_deep)

  nz = 400
  photoz_dist, scat[w_lrg].photoz_z, scat[w_lrg].photoz_zerr, $
    nz, zvals, pofz

  ;; Make some plots
  pclr = [!yellow,!red, !green, !blue, !magenta, !p.color,!seaGreen]
  lines = [0,0,0,0,2,0,0]
  names=['LRG low med','LRG hi med','LRG low deep','LRG high deep','All','All LRG','Lens Galaxies']
  xrange = [0,1]
  yrange = [0,6]
  xtitle = 'photoz'
  ytitle = 'P(photoz)'
  plot, zvals, pofz, yrange=yrange, xrange=xrange, xstyle=1+2,$
    xtitle=xtitle,ytitle=ytitle
  oplot,nzstruct.z,nzstruct.pofz,line=2,color=pclr[4]
  plothist,lcat.z,bin=0.01,norm=0.5,max=0.4,/overplot,color=pclr[6]

  legend, names,lines=lines,color=pclr,/right,box=0,charsize=1

  photoz_dist, $
    scat[w_lrg_lo_med].photoz_z, scat[w_lrg_lo_med].photoz_zerr, $
    nz, zvals, pofz
  fac = float(n_lrg_lo_med)/n_lrg
  oplot, zvals, pofz*fac, color=pclr[0]

  photoz_dist, $
    scat[w_lrg_hi_med].photoz_z, scat[w_lrg_hi_med].photoz_zerr, $
    nz, zvals, pofz
  fac = float(n_lrg_hi_med)/n_lrg
  oplot, zvals, pofz*fac, color=pclr[1]

  photoz_dist, $
    scat[w_lrg_lo_deep].photoz_z, scat[w_lrg_lo_deep].photoz_zerr, $
    nz, zvals, pofz
  fac = float(n_lrg_lo_deep)/n_lrg
  oplot, zvals, pofz*fac, color=pclr[2]

  photoz_dist, $
    scat[w_lrg_hi_deep].photoz_z, scat[w_lrg_hi_deep].photoz_zerr, $
    nz, zvals, pofz
  fac = float(n_lrg_hi_deep)/n_lrg
  oplot, zvals, pofz*fac, color=pclr[3]

return
END 
