PRO run_wtheta_lumw, wclr, ndup, stripe, indir=indir, outdir=outdir, $
                     use_lambda=use_lambda,$
                     numfile=numfile,nbeta=nbeta

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: run_wtheta_lumw, wclr, ndup, stripe, indir=indir,outdir=outdir, use_lambda=use_lambda'
      return
  ENDIF 

  ;cosmo=2 for lambda

  colors=['u','g','r','i','z']

  sdssidl_setup,/silent
  setup_mystuff
  IF n_elements(outdir) EQ 0 THEN $
    outdir = '/sdss5/data0/lensout/stripe'+ntostr(long(stripe))+'/'

  ;; this catalog contains "main" sample, with and
  ;; without spectral info
  get_tsgals,stripe,allcat,/spec,count=nall

  ;; z1d default is -1, so make cut to find 
  ;; matched spec gals
  ws=where( (allcat.z1d GT 0.0) ,ngal)
  spectra = allcat[ws]

  allcat=0

  CASE wclr OF
      0: file='/sdss5/data0/wtheta/wtheta-752-756-gal-allmag-u21.5.fit'
      1: file='/sdss5/data0/wtheta/wtheta-752-756-gal-allmag-g21.fit'
      2: file='/sdss5/data0/wtheta/wtheta-752-756-gal-allmag-r21.fit'
      3: file='/sdss5/data0/wtheta/wtheta-752-756-gal-allmag-i21.fit'
      4: file='/sdss5/data0/wtheta/wtheta-752-756-gal-allmag-z21.fit'
  ENDCASE 

  allcat = mrdfits(file, 1)

  ;; mag limit is to avoid overweighting bright galaxies
  wa=where(allcat.petrocounts[2] GT 15.5,nall)
  allcat = temporary(allcat[wa])
  IF stripe GT 45 THEN BEGIN 
      lambda=allcat.lambda
      rotate_lambda, lambda
      s=sort(lambda)
      lambda=0
  ENDIF ELSE s=sort(allcat.lambda)

  allcat = temporary(allcat[s])
;  nall = n_elements(allcat)

  print
  print,'Spectra Galaxies: ',ngal
  print
  print
  print,'All Galaxies: ',nall
  print

  rmin = 20.
  rmax = 1000.
  binsize = 40.

  print,'--------------------------'
  print,'Using ndup = ',ndup
  print,'--------------------------'
  print

  step = long(300*(ngal/30000.))

  wtheta_lensgal_lumweight, wclr, stripe, spectra, allcat, rmin, rmax, binsize, $
    use_lambda=use_lambda, step=step, $
    sumfile=sumfile, lensumfile=lensumfile, zfile=zfile, $
    outdir=outdir,$
    numfile=numfile,nbeta=nbeta

  tmp=mrdfits(zfile,1)
  zrand = tmp.z
  inputlambda = tmp.lambda

  numrand = ndup*n_elements(zrand)
  step = long(300*(numrand/30000.))
  
  wthetarand_lensgal_lumweight,wclr,stripe,allcat,rmin,rmax,binsize,zrand,ndup, $
    use_lambda=use_lambda, step=step, $
    outdir=outdir, zfile=zfile,$
    numfile=numfile,nbeta=nbeta

  delvarx, scat, tmp, zrand
  
  return
END 
