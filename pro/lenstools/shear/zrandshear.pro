
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    ZRANDSHEAR_LAMBDA
;       
; PURPOSE:
;    Measure the shear and density contrast for a set of lenses. Structure
;    must have a redshift tags.
;
; CALLING SEQUENCE:
;    zrandshear_lambda, stripe, lenscat, scat, clr, rmin, rmax, binsize, zrand, ndup, $
;                      use_lambda=use_lambda, $
;                      step=step, addstr=addstr, $
;                      outdir=outdir, $
;                      maxe=maxe, $
;                      wgood=wgood, $
;                      usecat=usecat, $
;                      datfile=datfile, sumfile=sumfile, zfile=zfile, $
;                      lensumfile=lensumfile
;
; INPUTS: 
;    stripe: stripe number in integer form.
;    lenscat: lens catalog
;    scat: source catalog. Must contain lambda and eta, and must be sorted
;          by lambda.
;    clr: the sdss bandpass
;    rmin: minimum radius in kpc
;    rmax: maximum radius in kpc
;    binsize: size of bins in kpc
;    zrand: the redshifts of the lenses spit out by zobjshear_lambda
;    ndup: now many times to duplicate each redshift.
;
; OPTIONAL INPUTS:
;    use_lambda: cosmology. omegamatter=1 is cosmo=1 (default). cosmo>1 is lambda=0.7
;    step: how many lenses to use in a group. default is 300
;    addstr: string to add to front of output file names.
;    outdir: output directory
;    maxe: maximum allowed ellipticity of source galaxy distribution.
;
; KEYWORD PARAMETERS:
;       
; OUTPUTS: 
;    datafile, sumfile, zfile, lensumfile
;
; OPTIONAL OUTPUTS:
;    wgood: the indices of used lenses.
;    usecat: the final lenses.
;    datafile, sumfile, zfile, lensumfile: the names of output files
;
; CALLED ROUTINES:
;    NTOSTR
;    TAG_EXIST
;    ZSHSTRUCT
;    ZSUMSTRUCT
;    ZLENSUMSTRUCT
;    EXIST
;    NEWNAME
;    (EQ2SURVEY)
;    (ROTATE_LAMBDA)
;    ARRSCL
;    ANGDIST_LAMBDA
;    SDSS_SIGMA_CRIT
;    ZRANDSHEAR_LAMBDA_GET_SOURCES (in this file)
;    BINARY_SEARCH
;    ZSHHDR
;    ZSUMHDR
;    MWRFITS
;    
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    03-DEC-2000: First documentation. Erin Scott Sheldon.
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO zrandshear_lambda_generate_lambdas, nlam, angmax, lambda

  COMMON rand_block, seed, FIRSTLAMBDA, LASTLAMBDA

  lambda = randomu(seed,nlam)

  FOR i=0L, nlam-1 DO BEGIN 

      lambda[i] = arrscl(lambda[i], $
                         FIRSTLAMBDA + angmax[i], $
                         LASTLAMBDA  - angmax[i], $
                         arrmin=0., arrmax=1.)
  ENDFOR 

END 


PRO zrandshear, stripe, scat, clr, rmin, rmax, $
                nbin_or_binsize, zrand_in, ndup, $
                lensra, lensdec, $
                $
                use_lambda=use_lambda, $
                step=step, addstr=addstr, $
                outdir=outdir, $
                wgood=wgood, $
                check=check, $
                datfile=datfile, sumfile=sumfile, zfile=zfile, $
                lensumfile=lensumfile, $
                maxe=maxe, maxit=maxit, $
                maskdir=maskdir, etarangedir=etarangedir, $
                commonsrc=commonsrc, logbin=logbin, southrot=southrot
  
  IF n_params() LT 10 THEN BEGIN
      print,'-Syntax: zrandshear, stripe, scat, clr, rmin, rmax, nbin_or_binsize, zrand, ndup,'
      print,'  use_lambda=use_lambda,'
      print,'  step=step, addstr=addstr, '
      print,'  outdir=outdir, '
      print,'  wgood=wgood, '
      print,'  check=check, '
      print,'  datfile=datfile, sumfile=sumfile, zfile=zfile,'
      print,'  lensumfile=lensumfile, '
      print,'  maxe=maxe, maxit=maxit, inputra=inputra, '
      print,'  maskdir=maskdir, etarangedir=etarangedir,'
      print,'  commonsrc=commonsrc'
      return
  ENDIF 

  colors = ['u','g','r','i','z']

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  COMMON rand_block, seed, FIRSTRA, LASTRA

  time = systime(1)

  IF stripe GT 45 THEN BEGIN 
      print
      print,'This is a Southern Stripe.'
      IF keyword_set(southrot) THEN BEGIN 
          print,'But it was rotated'
          issouth=0
      ENDIF ELSE issouth=1
      print
  ENDIF ELSE issouth=0

  vint = .32^2
  IF n_elements(maxe) EQ 0 THEN maxe = .2

  ;; Size by which to group things 300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 300L ELSE step = long(step)
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = ''
  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss4/data1/esheldon/TMP/'
  IF NOT keyword_set(check) THEN check=0
  
  stripestr = '_stripe'+ntostr(stripe)
  clr_str = '_'+colors[clr]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; maxit is number of retries until random point
  ;; passes symmetry check
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(maxit) EQ 0 THEN maxit = 100
                  

  d2r = !dpi/180.0d0            ;Change degrees to radians.
  r2d = 180./!dpi               ;radians to degrees

  ;;;;;;;;;;;;;;;;;;;;;;
  ;; log binning?
  ;;;;;;;;;;;;;;;;;;;;;;
  
  IF keyword_set(logbin) THEN BEGIN 
      print,'--------------------------'
      print,'Using Logarithmic Binning'
      print,'--------------------------'
      binsize = -1.
      nbin = nbin_or_binsize
  ENDIF ELSE BEGIN 
      binsize = nbin_or_binsize
      nbin = long( (rmax - rmin)/binsize ) ;+ 1 
  ENDELSE 

  print
  print,'-------------------------------------------------------------'
  print,'Rmin: ',rmin
  print,'Rmax: ',rmax
  print,'Step: ',step
  print,'Using ',ntostr(nbin),' bins between ',ntostr(rmin), $
        ' and ',ntostr(rmax),' kpc'

  IF NOT tag_exist(scat,'momerr',index=momerr) THEN BEGIN
      IF NOT tag_exist(scat,'uncert',index=momerr) THEN BEGIN 
          print
          print,'No valid moment uncertainty tag'
          print
          return
      ENDIF 
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; duplicate the input z distribution
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Duplicating: ',ntostr(ndup)

  nrand_in = n_elements(zrand_in)
  nrand = nrand_in*ndup

  ;; check input arrays
  nlensra = n_elements(lensra)
  nlensdec = n_elements(lensdec)

  IF (nlensra NE nrand_in) OR (nlensdec NE nrand_in) THEN BEGIN 
      message,'# of input ra not equal to # of input redshifts' 
  ENDIF ELSE BEGIN
      print & print,'Using input ra array' & print
  ENDELSE 

  FOR dd=0L, ndup-1 DO BEGIN 

      add_arrval, zrand_in, zrand
      add_arrval, lindgen(nrand_in), zindex

      ;; randomize ra/dec by generating random number
      ;; and sorting-> resulting indices are random
      r_ra=randomu(seed, nrand_in)
      r_dec=randomu(seed, nrand_in)

      s_ra=sort(r_ra)
      s_dec=sort(r_dec)

      add_arrval, lensra[s_ra], lra
      add_arrval, lensdec[s_dec], ldec

  ENDFOR 

  ;; randomize the redshifts
  tmprand = randomu(seed, nrand)
  s=sort(tmprand)
  zrand = zrand[s]
  zindex = zindex[s]

  IF issouth THEN rotate_ra, lra
  s=sort(lra)
  lra=lra[s]
  ldec=ldec[s]
  IF issouth THEN rotate_ra, lra

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; declare some arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; define structs and arrays

  lensumstruct = zlensumstruct(fltarr(nbin))
  lensumstruct = create_struct('z', 0., 'ra',0d, 'dec',0d, lensumstruct)

  lensum = replicate(lensumstruct, nrand )
  lensum.zindex = zindex
  lensum.z = zrand

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF addstr EQ '' THEN prename = outdir+addstr+'zrand'+stripestr+clr_str $
  ELSE prename = outdir+addstr+'_zrand'+stripestr+clr_str
  zobjshear_names, prename, datfile, sumfile, zfile, lensumfile, groupfile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; source cat
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; they are already sorted
  nsource = n_elements(scat)
  wsource = lindgen(nsource)

  FIRSTRA = scat[0].ra
  LASTRA  = scat[nsource-1].ra

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Lens cat: Set up sigma crit and Dlens. 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  h=1.
  print,'Using h = ',h
  print,'Calculating 1/sigma_crit'
  print

  wlens = lindgen(nrand)

  sigcritinv = sdss_sigma_crit(stripe, clr, zrand, $
                               wgood=wgood, use_lambda=use_lambda, $
                               commonsrc=commonsrc)
  sigmacrit = 1./sigcritinv
  wlens = wlens[wgood]
  lensum.scritinv = sigcritinv

  nlens = n_elements(wlens)

  IF NOT keyword_set(use_lambda) THEN omegamat=1.0 ELSE omegamat=0.3
  DL = angdist_lambda( zrand, h=h, omegamat=omegamat)*1000. ;Convert from Mpc to kpc
  angmax = rmax/DL*180./!pi     ;angmax in degrees

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Throw out by _edge_ if requested (rather than just symmety of background dist.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF check THEN BEGIN           ;May just want to find the good lenses.
      wgood = wlens
      return
  ENDIF 

  nlens = n_elements(wlens)
  print,'Threw out ',ntostr(nlens-nrand),' edge lenses'
  print
  print,'Using '+ntostr(nlens)+'/'+ntostr(nrand)+' lenses'
  print

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Estimate shear around all of these lenses
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  zobjshear_looplens, lensum, scat, lra, ldec, angmax, DL,$
    step, maxe, wlens, $
    rmin, rmax, nbin_or_binsize, $
    groupstruct, indices, $
    logbin=logbin, issouth=issouth, maxit=maxit,$
    /random

  ;; Remove unused lenses
  wbad = where(indices EQ -1, nwbad)
  IF nwbad NE 0 THEN remove, wbad, wlens
  lensused = n_elements(wlens)
  
  wgood = wlens

  print
  print,'Finally used ',ntostr(lensused),'/',ntostr(nlens),' Lenses'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find averages from sums
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  
  lensum = temporary(lensum[wlens])
  lensw = lensum.scritinv^2*lensum.totpairs
  lenswsum = total( lensw )
  zsum = total(zrand[wlens]*lensw)
  zmean = zsum/lenswsum
  print
  print,'zmean = ',zmean

  scritinvsum = total( lensum.scritinv*lensw )
  meanscritinv = scritinvsum/lenswsum
  print
  print,'meanscritinv: ',meanscritinv
  print

  combine_zlensum, lensum, binsize, rmin, rmax, h, sumstruct, shstruct

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Ouput the data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Shear for these lenses
  shhdr = zshhdr(shstruct)   
  mwrfits2, shstruct, datfile, shhdr, /create, /destroy

  ;; Sum file for these lenses
  sumhdr = zsumhdr(sumstruct)
  lhdr = sumhdr
  mwrfits2, sumstruct, sumfile, sumhdr, /create, /destroy

  ;; File with structure for each lens used
  lensum.ra = lra[wlens]
  lensum.dec = ldec[wlens]

  ;; File containing ra,dec,redshift for each used lens
  zs = create_struct('z', 0., $
                     'scritinv', 0., $
                     'ra', double(0.), $
                     'dec', double(0.) )
  zstruct = replicate(zs, lensused)
  zstruct.z = zrand[wlens]
  zstruct.scritinv = lensum.scritinv
  zstruct.ra = lra[wlens]
  zstruct.dec = ldec[wlens]

  mwrfits2, zstruct, zfile, /create, /destroy

  ;; file containing group stats
  mwrfits2, groupstruct, groupfile, /create, /destroy

  mwrfits2, lensum, lensumfile, lhdr, /create, /destroy

  step=oldstep

  ptime, systime(1)-time
  return
END 
