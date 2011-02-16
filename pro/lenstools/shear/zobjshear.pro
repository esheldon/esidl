
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    ZOBJSHEAR_LAMBDA
;       
; PURPOSE:
;    Measure the shear and density contrast for a set of lenses. Structure
;    must have a redshift tags.
;
; CALLING SEQUENCE:
;    zobjshear_lambda, stripe, lenscat, scat, clr, rmin, rmax, binsize, $
;                      use_lambda=use_lambda, $
;                      step=step, addstr=addstr, $
;                      outdir=outdir, $
;                      maxe=maxe, $
;                      wgood=wgood, $
;                      check=check, $
;                      datfile=datfile, sumfile=sumfile, zfile=zfile, $
;                      lensumfile=lensumfile
;
; INPUTS: 
;    stripe: stripe number in integer form.
;    lenscat: lens catalog. Must have either (ra,dec) or (lambda,eta)
;    scat: source catalog. Must contain (lambda,eta), and must be sorted
;          by lambda.
;    clr: the sdss bandpass
;    rmin: minimum radius in kpc
;    rmax: maximum radius in kpc
;    binsize: size of bins in kpc
;
; OPTIONAL INPUTS:
;    /use_lambda: lambda=0.7 else lambda=0 (omega_tot=1)
;    step: how many lenses to use in a group. default is 300
;    addstr: string to add to front of output file names.
;    outdir: output directory
;    maxe: maximum allowed ellipticity of source galaxy distribution.

; KEYWORD PARAMETERS:
;    /check: just see which lenses are in allowd redshift range, set wgood and return.
;       
; OUTPUTS: 
;    datafile, sumfile, zfile, lensumfile
;
; OPTIONAL OUTPUTS:
;    wgood: the indices of used lenses.
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
;    ANGDIST_LAMBDA
;    COPY_STRUCT
;    SDSS_SIGMA_CRIT
;    ZOBJSHEAR_LAMBDA_GET_SOURCES (in this file)
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


PRO zobjshear, stripe, lenscat, scat, clr, rmin, rmax, $
                      nbin_or_binsize, $
                      commonsrc=commonsrc,$
                      use_lambda=use_lambda, $
                      step=step, addstr=addstr, $
                      outdir=outdir, $
                      wgood=wgood, $
                      check=check, $
                      datfile=datfile, sumfile=sumfile, zfile=zfile, $
                      lensumfile=lensumfile, $
                      maxe=maxe, logbin=logbin, southrot=southrot

  IF n_params() LT 7 THEN BEGIN
      print,'-Syntax: zobjshear, stripe, lenscat, scat, clr, rmin, rmax, nbin,'
      print,'  use_lambda=use_lambda,'
      print,'  step=step, addstr=addstr, '
      print,'  outdir=outdir, '
      print,'  wgood=wgood, '
      print,'  check=check, '
      print,'  datfile=datfile, sumfile=sumfile, zfile=zfile, '
      print,'  lensumfile=lensumfile, '
      print,'  maxe=maxe '
      print,'-> lenscat may have ra,dec or lambda, eta but scat must contain lambda,eta'
      print,'-> scat must be sorted by lambda'
      return
  ENDIF 

  colors = ['u','g','r','i','z']

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

  ;; Size by which to group things  300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 300L
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = 'zgal_gal'
  IF n_elements(outdir) EQ 0 THEN outdir = '/sdss4/data1/esheldon/TMP/'
  IF NOT keyword_set(check) THEN check=0

  stripestr = '_stripe'+ntostr(stripe)
  clr_str = '_'+colors[clr]

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
  IF NOT keyword_set(use_lambda) THEN print,'Using lambda = 0' $
  ELSE print,'Using lambda = 0.3'
  print,'Rmin: ',rmin
  print,'Rmax: ',rmax
  print,'Step: ',step
  print,'Using ',ntostr(nbin),' bins between ',ntostr(rmin), $
        ' and ',ntostr(rmax),' kpc'

  ;; hubble parameter
  h=1.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; declare some arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; define structs and arrays. 

  lensumstruct = zlensumstruct(fltarr(nbin))
  lensumstruct = create_struct(lenscat[0], lensumstruct)

  lensum = replicate(lensumstruct, n_elements(lenscat) )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  prename = outdir+addstr+stripestr+clr_str
  zobjshear_names, prename, datfile, sumfile, zfile, lensumfile, groupfile,$
    indfile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; source cat
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; uncertainty tag
  momerr = geterrtag(scat[0])

  ;; they are already sorted
  nsource = n_elements(scat)
  wsource = lindgen(nsource)

  FIRSTRA = scat[0].ra
  LASTRA  = scat[nsource-1].ra

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Lens cat set up
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  copy_struct, lenscat, lensum

  zobjshear_setuplens, lensum, rmax, FIRSTRA, LASTRA, stripe, clr, $
    wz, lra, ldec, wlens, DL, angmax, $
    tsgals=tsgals, maskdir=maskdir, $
    use_lambda=use_lambda, commonsrc=commonsrc, h=h, issouth=issouth

  nlens = n_elements(wlens)

  ;; can return now if just checking
  IF check THEN BEGIN           ;May just want to find the good lenses.
      wgood = wlens
      return
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; print some stuff
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,'Using '+ntostr(nlens)+'/'+ntostr(n_elements(lensum))+' lenses'
  print
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Estimate shear around all of these lenses
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  zobjshear_looplens, lensum, scat, lra, ldec, angmax, DL,$
    step, maxe, wlens, $
    rmin, rmax, nbin_or_binsize, $
    groupstruct, indices, $
    logbin=logbin, issouth=issouth, maxit=maxit

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
  zsum = total(lensum.(wz[0])*lensw)
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

  ;; Shear file for these lenses
  shhdr = zshhdr(shstruct)      
  mwrfits2, shstruct, datfile, shhdr, /create, /destroy

  ;; Sum file for these lenses
  sumhdr = zsumhdr(sumstruct)
  lhdr = sumhdr
  mwrfits2, sumstruct, sumfile, sumhdr, /create, /destroy

  ;; File containing ra,dec,redshift for each used lens
  zs = create_struct('z', 0., $
                     'scritinv', 0., $
                     'ra', double(0.), $
                     'dec', double(0.) )
  zstruct = replicate(zs, lensused)
  zstruct.z = lensum.(wz[0])
  zstruct.scritinv = lensum.scritinv
  zstruct.ra = lra[wlens]
  zstruct.dec = ldec[wlens]

  mwrfits2, zstruct, zfile, /create, /destroy

  ;; file containing group stats
  mwrfits2, groupstruct, groupfile, /create, /destroy

  ;; File with structure for each lens used
;  lensum = temporary(lensum[wlens])
  lensum.zindex = lindgen(lensused)
  mwrfits2, lensum, lensumfile, lhdr, /create, /destroy

  step=oldstep

  ptime, systime(1)-time
  return
END 
