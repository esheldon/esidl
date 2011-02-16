PRO run_wtheta, ndup, stripe, indir=indir, outdir=outdir, use_lambda=use_lambda

  IF n_params() LT 2 THEN BEGIN
      print,'-Syntax: run_wtheta, ndup, stripe, indir=indir,outdir=outdir, use_lambda=use_lambda'
      return
  ENDIF 

  ;cosmo=2 for lambda

  colors=['u','g','r','i','z']

  sdssidl_setup,/silent
  setup_mystuff
  IF n_elements(outdir) EQ 0 THEN $
    outdir = '/sdss5/data0/lensout/stripe'+ntostr(long(stripe))+'/'

  ;; spectra lcat
;  get_spectra_lcat, stripe, spectra, indir=indir,count=ngal

  ;; "complete" catalog
;  get_lcat, stripe, 2, allcat, indir=indir,count=nall

  ;; this catalog contains "main" sample, with and
  ;; without spectral info
  get_tsgals,stripe,allcat,/spec,count=nall

  ;; z1d default is -1, so make cut to find 
  ;; matched spec gals
  ws=where(allcat.z1d GT 0.0,ngal)
  spectra = allcat[ws]

  eq2survey,allcat.ra,allcat.dec,lam,eta
  add_tags,allcat,['lambda','eta'],['0d','0d'],newstr
  newstr.lambda = lam
  newstr.eta = eta
  setzero,allcat
  IF stripe GT 45 THEN BEGIN 
      rotate_lambda,lam
      s=sort(lam)
  ENDIF ELSE s=sort(lam)
  allcat = temporary(newstr[s])

  print
  print,'Spectra Galaxies: ',ngal
  print
  print
  print,'All Galaxies: ',nall
  print

  rmin = 20.
  rmax = 1000.
  binsize = 80.

  print,'--------------------------'
  print,'Using ndup = ',ndup
  print,'--------------------------'
  print

  step = long(300*(ngal/30000.))
  wtheta_lensgal, stripe, spectra, allcat, rmin, rmax, binsize, $
    use_lambda=use_lambda, step=step, $
    sumfile=sumfile, lensumfile=lensumfile, zfile=zfile, $
    outdir=outdir

  tmp=mrdfits(zfile,1)
  zrand = tmp.z
  inputlambda = tmp.lambda

  numrand = ndup*n_elements(zrand)
  step = long(300*(numrand/30000.))
  
  wthetarand_lensgal, stripe, allcat, rmin, rmax, binsize, zrand, ndup, $
    use_lambda=use_lambda, step=step, $
    outdir=outdir, zfile=zfile

  delvarx, scat, tmp, zrand
  
  return
END 
