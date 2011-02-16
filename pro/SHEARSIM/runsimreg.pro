PRO runsimreg, eend, lcat=lcat, scat=scat, sigma=sigma, cutoff=cutoff,$
               use=use, nocut=nocut, snoise=snoise

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: runsimreg, sigma, cutoff, lcat, scat, eend, use=use, nocut=nocut'
      return
  ENDIF 

  IF n_elements(cutoff) EQ 0 THEN cutoff = 1200.
  IF n_elements(sigma) EQ 0 THEN sigma = 170.

;  zsource = .4
;  zlens = .15
  zsource = .285                ;makes sigcrit right...
  zlens = .172
  simgal_gal, sigma, cutoff, zsource, zlens, scat, lcat, $
    LASTRA=LASTRA, FIRSTRA=FIRSTRA, use=use, nocut=nocut, snoise=snoise

;  scat=0
;  sra=1.
;  sdec=1.
;  sis_shear, zlens, zsource, sigma, sdec, sra, scat, cen, h=1., omega=1.

;  lcat = replicate(lcat[0],1)
;  lcat.dec = cen[0]
;  lcat.ra = cen[1]

  run1=11
  run2=22
  clr=2
  rmin = 10.
  rmax = 600.
  binsize = 39.

  nbin = 16

  outdir = '/sdss4/data1/esheldon/MANYSIM/'
  addstr='sim'

  objshear_fix, run1, run2, clr, rmin, rmax, binsize, $
    FIRSTRA=FIRSTRA, LASTRA=LASTRA, scat=scat, lcat=lcat, $
    addstr = addstr, outdir=outdir

;  objshear, run1, run2, clr, rmin, rmax, binsize, $
;    scat=scat, lcat=lcat, outdir=outdir

  regressgal, run1, run2, clr, rmin, rmax, binsize, $
    FIRSTRA=FIRSTRA, LASTRA=LASTRA, scat=scat, lcat=lcat, $
    addstr = addstr, outdir=outdir

  matf=outdir+'simmatrix_11_22_r_N'
  evecf=outdir+'simevec_11_22_r_N'
  rf=outdir+'simrsum_11_22_r_N'
  
  IF n_elements(eend) EQ 0 THEN BEGIN 
      eend=' '
      print,format='($,"Enter N")'
      read,eend
  ENDIF 

  matf=matf+eend+'.txt'
  evecf=evecf+eend+'.txt'
  rf=rf+eend+'.txt'

  wm_trans_m, matf, evecf, ata, btb, nbin
  wsolve_reg,matf,evecf,nbin,ata,btb,etan1,erad1,etan2,erad2

  shear1 = etan1/2.
  ortho1 = erad1/2.
  shear2 = etan2/2.
  ortho2 = erad2/2.

  uncert1 = fltarr(nbin)
  uncert2 = fltarr(nbin)

  a_inv=invert(ata)
  b_inv=invert(btb)

  FOR i=0, nbin-1 DO BEGIN 
      uncert1[i] = sqrt(a_inv[i,i])/2.
      uncert2[i] = sqrt(b_inv[i,i])/2.
  ENDFOR 

  readcol, rf, meanr

  tmp=str_sep(matf, 'simmatrix')

  outfile = tmp[0]+'simregshears'+tmp[1]
  atafile = tmp[0]+'ata'+tmp[1]
  btbfile = tmp[0]+'btb'+tmp[1]
  IF (exist(atafile) OR exist(btbfile) OR exist(outfile) )THEN BEGIN 
      atafile = newname(atafile)
      btbfile = newname(btbfile)
      outfile = newname(outfile)
  ENDIF 
  print,'ata file: ',atafile
  print,'btb file: ',btbfile
  print,'outfile: ',outfile

  colprint, meanr, shear1, uncert1, shear2, uncert2, ortho1, ortho2
  
  fmt='(7(F0,:,1X))'
  openw, unit, outfile, /get_lun
  colprint,meanr, shear1, uncert1, shear2, uncert2, ortho1, ortho2, $
           lun=unit, format=fmt
  free_lun, unit
  mwrfits, ata, atafile, /create
  mwrfits, btb, btbfile, /create

  return

  rdobjshear, outdir+'sim11_22_r_N'+eend+'.dat',st1
;  wait,3
;  rdobjshear, outdir+'gal_gal_11_22_r_N'+eend+'.dat',st2

  sh=shearsis_trunc(sigma,cutoff, zsource,zlens,st1.meanr,/core)

  plot, st1.meanr, st1.shear, psym=1
;  oplot, st2.meanr, st2.shear, psym=2
  oplot, st1.meanr, sh
  oplot, meanr, shear1, line=1
  oplot, meanr, shear2, line=2

return
END 
