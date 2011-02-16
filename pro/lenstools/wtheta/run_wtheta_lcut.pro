PRO run_wtheta_lcut, wclr, ndup, stripe, indir=indir, outdir=outdir, $
                     use_lambda=use_lambda,fraclstar=fraclstar,$
                     numfile=numfile,nbeta=nbeta, tsgals=tsgals,rmax=rmax, $
                     edgecut=edgecut, maskdir=maskdir, $
                     etarangedir=etarangedir, randomize=randomize,$
                     zfile=zfile

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: run_wtheta_lcut, wclr, ndup, stripe, indir=indir, outdir=outdir, use_lambda=use_lambda,numfile=numfile,nbeta=nbeta, rmax=rmax, tsgals=tsgals,rmax=rmax'
      return
  ENDIF 

  ;cosmo=2 for lambda

  colors=['u','g','r','i','z']

  sdssidl_setup,/silent
  setup_mystuff

  stripestr = ntostr(fix(stripe))
  IF n_elements(outdir) EQ 0 THEN $
    outdir = '/sdss5/data0/lensout/stripe'+stripestr+'/'

  ;; this catalog contains "main" sample, with and
  ;; without spectral info
;  get_tsgals,stripe,acat,/spec,count=nall
  get_spectra_lcat, stripe, acat, count=nall


  ;; z1d default is -1, so make cut to find 
  ;; matched spec gals
  make_tsflag_struct, ts
  ts.galaxy='Y'
  tsflag_select, acat, ts, si
;  ws=where( (acat[si].z1d GT 0.0) ,ngal)
  ws = where( (acat[si].z1d GT 0.05),ngal)
  spectra = acat[ si[ws] ]

  acat=0

  ;; get rid of duplicates.  Why are they there?
  phid=photoid(spectra.run,spectra.rerun,spectra.camcol,spectra.field,spectra.id)
  spectra = temporary( spectra[rem_dup(phid)] )

;  windir='/sdss5/data0/wtheta/
;  CASE wclr OF
;      0: file=windir+'wtheta-752-756-gal-allmag-u21.5.fit'
;      1: file=windir+'wtheta-752-756-gal-allmag-g21.fit'
;      2: file=windir+'wtheta-752-756-gal-allmag-r21.fit'
;      3: file=windir+'wtheta-752-756-gal-allmag-i21.fit'
;      4: file=windir+'wtheta-752-756-gal-allmag-z21.fit'
;  ENDCASE 

  windir='/sdss6/data0/wtheta/'
  CASE wclr OF
      0: file=windir+'wtheta-stripe'+stripestr+'-gal-allmag3-u21.fit'
      1: file=windir+'wtheta-stripe'+stripestr+'-gal-allmag3-g21.fit'
      2: file=windir+'wtheta-stripe'+stripestr+'-gal-allmag3-r21.fit'
      3: file=windir+'wtheta-stripe'+stripestr+'-gal-allmag3-i21.fit'
      4: file=windir+'wtheta-stripe'+stripestr+'-gal-allmag3-z21.fit'
  ENDCASE 

  allcat = mrdfits(file, 1)

  ;; give random positions
  ;; for stripe 82 will have to rotate_ra first, then back
  IF keyword_set(randomize) THEN BEGIN 
      randomize_radec,allcat
      eq2survey, allcat.ra, allcat.dec, lam,eta
      allcat.lambda=lam
      allcat.eta = eta
      ss=sort(allcat.lambda)
      allcat = temporary(allcat[ss])

      randomize_radec,spectra
  ENDIF 

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

  print
  print,'Spectra Main Galaxies: ',ngal
  print
  print
  print,'All Galaxies: ',nall
  print

  rmin = 20.
  IF n_elements(rmax) EQ 0 THEN rmax = 1000.
  binsize = 40.

  print,'--------------------------'
  print,'Using ndup = ',ndup
  print,'--------------------------'
  print

  IF n_elements(zfile) NE 0 THEN GOTO,jump


  step = long(300*(ngal/30000.))

  wtheta_lensgal_lcut, wclr, stripe, spectra, allcat, rmin, rmax, binsize, $
    use_lambda=use_lambda, step=step, $
    sumfile=sumfile, lensumfile=lensumfile, zfile=zfile, $
    outdir=outdir,$
    numfile=numfile,nbeta=nbeta,fraclstar=fraclstar, tsgals=tsgals, $
    edgecut=edgecut, maskdir=maskdir, etarangedir=etarangedir

  IF ndup EQ -1 THEN return

jump:

  tmp=mrdfits(zfile,1)
  zrand = tmp.z
  inputlambda = tmp.lambda

  numrand = ndup*n_elements(zrand)
  step = long(300*(numrand/30000.))
  
  wthetarand_lensgal_lcut,wclr,stripe,allcat,rmin,rmax,binsize,zrand,ndup, $
    use_lambda=use_lambda, step=step, $
    outdir=outdir, zfile=zfile,$
    numfile=numfile,nbeta=nbeta,fraclstar=fraclstar, tsgals=tsgals, $
    edgecut=edgecut, maskdir=maskdir, etarangedir=etarangedir

  delvarx, scat, tmp, zrand

  return
END 
