PRO corrlowdense

  indir = '/sdss4/data1/esheldon/GAL_GAL/DENSITY/'

  ext = '_N1.dat'

  Ssh = 0.8777

  simpctable, rct, gct, bct

  !p.background = !white
  !p.color = !black

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; dense stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  dgf = indir + 'dense_752_756_g_N1.dat'
  drf = indir + 'dense_752_756_r_N1.dat'
  dif = indir + 'dense_752_756_i_N1.dat'

  rdobjshear, dgf, dgstruct, dgnlens,/silent
  rdobjshear, drf, drstruct, drnlens,/silent
  rdobjshear, dif, distruct, dinlens,/silent

  dgrndf = indir + 'dense_rand_752_756_g_N1.dat'
  drrndf = indir + 'dense_rand_752_756_r_N1.dat'
  dirndf = indir + 'dense_rand_752_756_i_N1.dat'

  rdobjshear, dgrndf, dgrndstruct, dgrndnlens  ,/silent
  rdobjshear, drrndf, drrndstruct, drrndnlens,/silent
  rdobjshear, dirndf, dirndstruct, dirndnlens,/silent

  ;; Don't use last bin
  nd =  n_elements(dgstruct.meanr)
  dgstruct = dgstruct[0:nd-2]
  drstruct = drstruct[0:nd-2]
  distruct = distruct[0:nd-2]

  ;; find correction
  dgL = dgstruct.npair/dgnlens
  drL = drstruct.npair/drnlens
  diL = distruct.npair/dinlens
  
  dgrndL = dgrndstruct.npair/dgrndnlens
  drrndL = drrndstruct.npair/drrndnlens
  dirndL = dirndstruct.npair/dirndnlens

  dgfrac =(dgL - dgrndL)/dgrndL > 0.
  drfrac =(drL - drrndL)/drrndL > 0.
  difrac =(diL - dirndL)/dirndL > 0.

  dgcorr = 1. + dgfrac
  drcorr = 1. + drfrac
  dicorr = 1. + difrac

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; low dense stuff
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  lgf = indir + 'low_752_756_g_N1.dat'
  lrf = indir + 'low_752_756_r_N1.dat'
  lif = indir + 'low_752_756_i_N1.dat'

  rdobjshear, lgf, lgstruct, lgnlens  ,/silent
  rdobjshear, lrf, lrstruct, lrnlens,/silent
  rdobjshear, lif, listruct, linlens,/silent

  lgrndf = indir + 'low_rand_752_756_g_N1.dat'
  lrrndf = indir + 'low_rand_752_756_r_N1.dat'
  lirndf = indir + 'low_rand_752_756_i_N1.dat'

  rdobjshear, lgrndf, lgrndstruct, lgrndnlens  ,/silent
  rdobjshear, lrrndf, lrrndstruct, lrrndnlens,/silent
  rdobjshear, lirndf, lirndstruct, lirndnlens,/silent

  ;; Don't use last bin
  nd =  n_elements(lgstruct.meanr)
  lgstruct = lgstruct[0:nd-2]
  lrstruct = lrstruct[0:nd-2]
  listruct = listruct[0:nd-2]

  ;; find correction
  lgL = lgstruct.npair/lgnlens
  lrL = lrstruct.npair/lrnlens
  liL = listruct.npair/linlens
  
  lgrndL = lgrndstruct.npair/lgrndnlens
  lrrndL = lrrndstruct.npair/lrrndnlens
  lirndL = lirndstruct.npair/lirndnlens

  lgfrac =(lgL - lgrndL)/lgrndL > 0.
  lrfrac =(lrL - lrrndL)/lrrndL > 0.
  lifrac =(liL - lirndL)/lirndL > 0.

  lgcorr = 1. + lgfrac
  lrcorr = 1. + lrfrac
  licorr = 1. + lifrac

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Correct shear
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  dgstruct.shear = dgstruct.shear/Ssh*dgcorr
  drstruct.shear = drstruct.shear/Ssh*drcorr
  distruct.shear = distruct.shear/Ssh*dicorr

  lgstruct.shear = lgstruct.shear/Ssh*lgcorr
  lrstruct.shear = lrstruct.shear/Ssh*lrcorr
  listruct.shear = listruct.shear/Ssh*licorr

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output corrected shear
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; high density
  dgf_out = indir + 'shear_dense_corr_g.txt'
  drf_out = indir + 'shear_dense_corr_r.txt'
  dif_out = indir + 'shear_dense_corr_i.txt'

  openw, lun, dgf_out, /get_lun
  colprint, dgstruct.meanr, dgstruct.shear, dgstruct.shearerr, lun=lun
  free_lun, lun

  openw, lun, drf_out, /get_lun
  colprint, drstruct.meanr, drstruct.shear, drstruct.shearerr, lun=lun
  free_lun, lun

  openw, lun, dif_out, /get_lun
  colprint, distruct.meanr, distruct.shear, distruct.shearerr, lun=lun
  free_lun, lun

  ;; low density
  lgf_out = indir + 'shear_low_corr_g.txt'
  lrf_out = indir + 'shear_low_corr_r.txt'
  lif_out = indir + 'shear_low_corr_i.txt'

  openw, lun, lgf_out, /get_lun
  colprint, lgstruct.meanr, lgstruct.shear, lgstruct.shearerr, lun=lun
  free_lun, lun

  openw, lun, lrf_out, /get_lun
  colprint, lrstruct.meanr, lrstruct.shear, lrstruct.shearerr, lun=lun
  free_lun, lun

  openw, lun, lif_out, /get_lun
  colprint, listruct.meanr, listruct.shear, listruct.shearerr, lun=lun
  free_lun, lun

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; plot dense g,r,i
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;

  yrange=fltarr(2)
  drange = prange(drstruct.shear, distruct.shear, drstruct.shearerr,distruct.shearerr)
  lrange = prange(lrstruct.shear, listruct.shear, lrstruct.shearerr,listruct.shearerr)
  yrange[0] = min( [drange[0], lrange[0]] )
  yrange[1] = max( [drange[1], lrange[1]] )

  erase & multiplot, [1,3]

  ploterror, dgstruct.meanr, dgstruct.shear, dgstruct.shearerr, psym=1, yrange=yrange
  oplot,[0,10000],[0,0]

  multiplot
  ploterror, drstruct.meanr, drstruct.shear, drstruct.shearerr, psym=1, yrange=yrange
  oplot,[0,10000],[0,0]

  multiplot
  ploterror, distruct.meanr, distruct.shear, distruct.shearerr, psym=1, yrange=yrange, $
    xtitle='Projected Radius (arcsec)', ytitle='Tangential Shear'
  oplot,[0,10000],[0,0]

  multiplot,/reset

  key = get_kbrd(1)

  ;; plot low dense g,r,i

  erase & multiplot, [1,3]

  ploterror, lgstruct.meanr, lgstruct.shear, lgstruct.shearerr, psym=1, yrange=yrange
  oplot,[0,10000],[0,0]

  multiplot
  ploterror, lrstruct.meanr, lrstruct.shear, lrstruct.shearerr, psym=1, yrange=yrange
  oplot,[0,10000],[0,0]

  multiplot
  ploterror, listruct.meanr, listruct.shear, listruct.shearerr, psym=1, yrange=yrange,$
    xtitle='Projected Radius (arcsec)', ytitle='Tangential Shear'
  oplot,[0,10000],[0,0]

  multiplot,/reset

  key = get_kbrd(1)

  ;; Plot together

  erase & multiplot, [1,3]

  plot, [0], /nodata, yrange=yrange, xrange=[0,600]
  oploterror, lgstruct.meanr, lgstruct.shear, lgstruct.shearerr, line=2, $
    color=!blue, errcolor=!blue
  oploterror, dgstruct.meanr, dgstruct.shear, dgstruct.shearerr, line=0, $
    color=!red, errcolor=!red
  oplot,[0,10000],[0,0]

  multiplot
  plot, [0], /nodata, yrange=yrange, xrange=[0,600]
  oploterror, lrstruct.meanr, lrstruct.shear, lrstruct.shearerr, line=2, $
    color=!blue, errcolor=!blue
  oploterror, drstruct.meanr, drstruct.shear, drstruct.shearerr, line=0, $
    color=!red, errcolor=!red
  oplot,[0,10000],[0,0]

  multiplot
  plot, [0], /nodata, yrange=yrange, xrange=[0,600], $
    xtitle='Projected Radius (arcsec)', ytitle='Tangential Shear'
  oploterror, listruct.meanr, listruct.shear, listruct.shearerr, line=2, $
    color=!blue, errcolor=!blue
  oploterror, distruct.meanr, distruct.shear, distruct.shearerr, line=0, $
    color=!red, errcolor=!red
  oplot,[0,10000],[0,0]

  legend,['High Density', 'Low Density'],line=[0,2],colors=[!red,!blue],/right
  

  multiplot,/reset

  key=get_kbrd(1)

  ;; Just i-band

  !p.thick = 2
  !p.charthick = 2
  !p.charsize = 1.5
  !x.thick = 2
  !y.thick = 2

  ;; Note : always use r-band conv, its most complete!!!!!!
  innddir='/sdss4/data1/esheldon/TMP/conv/'
  ;; This was chosen after the fact from best-fit
  innffile=innddir+'conv_dense_r_cut600.txt'
  readcol, innffile, modelx, kappa, kappa_int, /silent

  ;; need fudge factor because I made my convolutions 
  ;; with these redshifts.  This corrects to right sigcrit
  ;; should do sigmacrit right.
  zs = .4
  zl = .15
  philzl = [0., 0.168, 0.172, 0.173, 0.]
  phil_sigcrit = 1./[1., 0.343, 0.392, 0.403, 1.]
  my_sigcrit = sigmacrit(zs, zl, h=1.)
  philDl = angdist(philzl, h=1.)
  myDl = angdist(zl, h=1.)
  fixfac = my_sigcrit*myDl/phil_sigcrit[2]/philDl[2]

  shear_sig = (kappa_int - kappa)*fixfac*(170.)^2

  plot, [0], /nodata, yrange=yrange, xrange=[0,600], $
    xtitle='Projected Radius (arcsec)', ytitle='Tangential Shear',$
    title='High/Low Density Regions'
  oploterror, listruct.meanr, listruct.shear, listruct.shearerr, line=2, $
    color=!blue, errcolor=!blue
  oploterror, distruct.meanr, distruct.shear, distruct.shearerr, line=0, $
    color=!red, errcolor=!red
  oplot,[0,10000],[0,0]
  oplot, modelx, shear_sig, color=black,line=3
  
  !p.charthick = 1
  legend,['High Density', 'Low Density','Neighbors'],$
    line=[0,2,3],colors=[!red,!blue,!black],/right,thick=[2,2,2]
  !p.charthick = 2

  write_gif, '/sdss4/data1/esheldon/NICEPLOTS/high_low_shear.gif',tvrd(),rct,gct,bct

  key = get_kbrd(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; contour plots for combined data
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ;; This is for combining all bands
  gg = 3.50959599
  rr = 2.235299689
  ii = 2.549197961
  gr = 1.586466533 
  gi = 1.6983151 
  ri = 1.512359722

  xtitle='Cutoff Radius (arcsec)'
  ytitle='!7r!3!DV!N (km/s)'

  ;; higher density
  ;; just use r-band model, most complete!!!!!
  clr = 2
  type = [0,1,2]                ;type 0-dense 1-low 2-tot
  rdwtheta_conv, clr, type[0], datax, data, dataerr, modelx, model, $
    sigma, cutoff

  title='Galaxy Parameters   High Density Regions'
  xrange = [min(cutoff),max(cutoff)]

  chisq_conf_3band, $
    dgstruct.meanr, dgstruct.shear, dgstruct.shearerr, $
    drstruct.meanr, drstruct.shear, drstruct.shearerr, $
    distruct.meanr, distruct.shear, distruct.shearerr, $
    gg, rr, ii, gr, gi, ri, $
    modelx, model, cutoff, sigma, $
    chisq_surf, min1, min2, low1, high1, low2, high2, $
    xtitle=xtitle, title=title,ytitle=ytitle, $
    xrange=xrange, xstyle=1, c_colors=!red

  write_gif, '/sdss4/data1/esheldon/NICEPLOTS/high_conf.gif',tvrd(),rct,gct,bct
  key = get_kbrd(1)

  ;; lower density
  ;; just use r-band model, most complete!!!!!
  clr = 2
  type = [0,1,2]                ;type 0-dense 1-low 2-tot
  rdwtheta_conv, clr, type[1], datax, data, dataerr, modelx, model, $
    sigma, cutoff

  title='Galaxy Parameters   Low Density Regions'
  xrange = [min(cutoff),max(cutoff)]

  chisq_conf_3band, $
    lgstruct.meanr, lgstruct.shear, lgstruct.shearerr, $
    lrstruct.meanr, lrstruct.shear, lrstruct.shearerr, $
    listruct.meanr, listruct.shear, listruct.shearerr, $
    gg, rr, ii, gr, gi, ri, $
    modelx, model, cutoff, sigma, $
    chisq_surf, min1, min2, low1, high1, low2, high2, $
    xtitle=xtitle, title=title,ytitle=ytitle, $
    xrange=xrange, xstyle=1, c_colors=!blue

  write_gif, '/sdss4/data1/esheldon/NICEPLOTS/low_conf.gif',tvrd(),rct,gct,bct

  !p.thick = 1
  !x.thick = 1
  !y.thick = 1
  !p.charsize = 1.
  !p.charthick=1



  return
END 
