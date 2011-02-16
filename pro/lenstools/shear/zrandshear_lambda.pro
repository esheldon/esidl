
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
;                      usecat=usecat, $
;                      datfile=datfile, sumfile=sumfile, zfile=zfile, $
;                      lensumfile=lensumfile, hirata=hirata
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

PRO zrandshear_lambda, stripe, scat, clr, rmin, rmax, $
                       nbin_or_binsize, zrand_in, ndup, $
                       $
                       use_lambda=use_lambda, $
                       step=step, addstr=addstr, $
                       outdir=outdir, $
                       check=check, $
                       datfile=datfile, sumfile=sumfile, zfile=zfile, $
                       lensumfile=lensumfile, $
                       maxe=maxe, maxit=maxit, $
                       logbin=logbin,$
                       depth=depth, photoz=photoz, recorr=recorr, $
                       extno=extno, compcut=compcut, combined=combined, $
                       hirata=hirata, print2ascii=print2ascii, $
                       clusters=clusters, lrg_sources=lrg_sources
  
  IF n_params() LT 8 THEN BEGIN
      print,'-Syntax: '
      return
  ENDIF 

  colors = ['u','g','r','i','z']

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time = systime(1)

  vint = .32^2
  IF n_elements(maxe) EQ 0 THEN maxe = .2

  IF n_elements(depth) EQ 0 THEN depth=10

  ;; Size by which to group things 300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 300L ELSE step = long(step)
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = ''
  IF n_elements(outdir) EQ 0 THEN $
    outdir = esheldon_config("lensout_dir")+'tmp/'
  IF NOT keyword_set(check) THEN check=0
  
  stripestr = '_stripe'+stripearr2string(stripe)
  clr_str = clrarr2string(clr)

  IF keyword_set(recorr) THEN recorrstr='_recorr' ELSE recorrstr=''
  IF keyword_set(hirata) THEN hiratastr = '_h' ELSE hiratastr=''
  IF keyword_set(lrg_sources) THEN lrgstr = '_lrg' ELSE lrgstr=''

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; maxit is number of retries until random point
  ;; passes symmetry check
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(maxit) EQ 0 THEN maxit = 100
                  
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

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; duplicate the input z distribution
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Duplicating: ',ntostr(ndup)

  nrand_in = n_elements(zrand_in)
  nrand = nrand_in*ndup

  FOR dd=0L, ndup-1 DO BEGIN 

      add_arrval, zrand_in, zrand
      add_arrval, lindgen(nrand_in), zindex

  ENDFOR 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; declare some arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; define structs and arrays

  lensumstruct = zlensumstruct(fltarr(nbin))
  lensumstruct = create_struct('z', 0., 'clambda',0d, 'ceta',0d, lensumstruct)

  lensum = replicate(lensumstruct, nrand )
  lensum.zindex = zindex
  lensum.z = zrand

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; if extno is sent, then does not check for existence of files
  IF addstr EQ '' THEN BEGIN
      prename = $
        outdir+addstr+ 'zrand'+stripestr+'_'+clr_str+lrgstr+recorrstr+hiratastr
  ENDIF ELSE BEGIN
      prename = $
        outdir+addstr+'_zrand'+stripestr+'_'+clr_str+lrgstr+recorrstr+hiratastr
  ENDELSE 

  IF keyword_set(print2ascii) THEN BEGIN 
      zobjshear_names, $
        prename, datfile, sumfile, zfile, lensumfile, groupfile,$
        indfile, psfile, asciifile, extno=extno
      openw, outlun, asciifile, /get_lun
  ENDIF ELSE BEGIN 
      zobjshear_names, $
        prename, datfile, sumfile, zfile, lensumfile, groupfile,$
        indfile, psfile, extno=extno
  ENDELSE 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; source cat
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  get_nz_photoz, stripe, clr, nzstruct, hirata=hirata, lrg_sources=lrg_sources
  wsource = nzstruct.useind
  IF keyword_set(recorr) THEN BEGIN 
      wsource2 = where( abs(scat[wsource].e1_recorr) LT 2 AND $
                        abs(scat[wsource].e2_recorr) LT 2 )
  ENDIF ELSE BEGIN 
      wsource2 = where( abs(scat[wsource].e1) LT 2 AND $
                        abs(scat[wsource].e2) LT 2 )
  ENDELSE 
  help,wsource,wsource2
  wsource = wsource[wsource2]

  ;; apply the seeing mask?
  IF keyword_set(combined) THEN BEGIN 
      clambda = scat[wsource].clambda
      ceta = scat[wsource].ceta

      print,' **> Making seeing and other cuts to source catalog'
      apply_pixel_mask,clambda,ceta,pbad,pgood,/combined

      IF pgood[0] NE -1 THEN wsource = wsource[pgood]
  ENDIF  

  nsource = n_elements(wsource)

  FIRSTLAMBDA = scat[wsource[0]].clambda
  LASTLAMBDA  = scat[wsource[nsource-1]].clambda 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; if /photoz, get scinv_struct to interpolate
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; this already includes the prior
  IF keyword_set(photoz) THEN BEGIN 
      get_meanscinv, stripe, clr, scinv_struct, use_lambda=use_lambda,$
                     hirata=hirata, lrg_sources=lrg_sources
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Lens cat: Set up sigma crit and Dlens. 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  hubble=1.
  print,'Using h = ',hubble
  print,'Calculating 1/sigma_crit'
  print

  sigcritinv = sdss_sigma_crit(stripe, clr, zrand, $
                               wgood=wgood, use_lambda=use_lambda,hirata=hirata)
  IF n_elements(wgood) NE nrand THEN message,'Invalid input redshifts'
  sigmacrit = 1./sigcritinv
  lensum.scritinv = sigcritinv

  max_allow_angle = 6.0         ;degrees
  IF NOT keyword_set(use_lambda) THEN omegamat=1.0 ELSE omegamat=0.3
  ;; angular diameter distance: Convert from Mpc to kpc
  DL = angdist_lambda( zrand, h=hubble, omegamat=omegamat)*1000. 
  angmax = rmax/DL*180.0/!pi     ;; angmax in degrees

  wbang = where(angmax GT max_allow_angle, nbang)
  IF nbang NE 0 THEN message,'Invalid redshifts: angle > 6.0 degrees'

  print,'Using '+ntostr(nrand)+' random lenses'
  print

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; plot out the redshifts and positions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  begplot,name=psfile

  !p.multi=[0,0,2]

  plothist, zrand, bin=0.01, xtitle='Z!DLens!N', ytitle='dN/dZ'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Estimate shear. Generate random points as needed
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  bad  = lindgen(nrand)
  nbad  = nrand
  ngood = 0
  good = 0L

  llambda = dblarr(nrand)
  leta    = llambda

  scat = temporary(scat[wsource])
  ;; Keep looping until we get nrand

  IF NOT keyword_set(clusters) THEN specgal = 1

  WHILE ngood LT nrand DO BEGIN 

      sdss_genrand, stripe, nbad, tlambda, teta, $
                    specgal=specgal, clusters=clusters, $
                    compcut=compcut, $
                    /photmask, angmax=angmax[bad], /twoquad, /nomissf, $
                    maskflags=maskflags, completeness=completeness
                    
      
      ;; flush stdout,stderr (doesn't help)
      ;;flush,-1,-2

      ;; copy in values (may change later)
      llambda[bad] = tlambda
      leta[bad]    = teta
      lensum[bad].clambda = tlambda
      lensum[bad].ceta = teta
      lensum[bad].pixelmaskflags = maskflags
      IF NOT keyword_set(clusters) THEN $
        lensum[bad].completeness   = completeness
      
      ;; these are sums over bins.  If already been through the code they
      ;; will have non-zero values and be added to
      lensum[bad].totpairs = 0
      lensum[bad].sshsum = 0
      lensum[bad].wsum_ssh = 0
      lensum[bad].weight = 0

      ;; note: sending bad
      zobjshear_htm_looplens, depth, lensum, scat, llambda, leta, angmax, DL,$
                              step, maxe, bad, $
                              rmin, rmax, nbin_or_binsize, $
                              groupstruct, indices, $
                              logbin=logbin, $
                              maxit=maxit,$
                              recorr=recorr,$
                              scinv_struct=scinv_struct, $
                              outlun=outlun

      ;; Remove unused lenses
      tbad = where(indices EQ -1, nbad, comp=tgood, ncomp=ntgood)
      
      IF ntgood NE 0 THEN BEGIN 
          ngood = ngood + ntgood
          IF nbad NE 0 THEN bad = bad[tbad] ELSE bad = -1L
      ENDIF 

  ENDWHILE 

  IF keyword_set(print2ascii) THEN free_lun, outlun

;  lensum.clambda = llambda
;  lensum.ceta    = leta

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find averages from sums
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  lensw = lensum.weight
  lenswsum = total( lensw )
  zsum = total(zrand*lensw)
  zmean = zsum/lenswsum
  print
  print,'zmean = ',zmean

  combine_zlensum, lensum, binsize, rmin, rmax, hubble, shstruct, compcut=compcut

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Finish plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  plot, llambda, leta, psym=3, /ynozero, $
        xtitle = !csym.lambda+'!Dc!N', ytitle = !csym.eta+'!Dc!N'
  !p.multi=0
  endplot 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Ouput the data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Shear for these lenses
  shhdr = zshhdr(shstruct)  
  lhdr  = zsumhdr(shstruct)
  mwrfits2, shstruct, datfile, shhdr, /create, /destroy

  ;; Sum file for these lenses
;  sumhdr = zsumhdr(sumstruct)
;  lhdr = sumhdr
;  mwrfits2, sumstruct, sumfile, sumhdr, /create, /destroy

  ;; File containing ra,dec,redshift for each used lens
  zs = create_struct('z', 0., $
                     'scritinv', 0., $
                     'clambda', double(0.), $
                     'ceta', double(0.) )
  zstruct = replicate(zs, nrand)
  zstruct.z = zrand
  zstruct.scritinv = lensum.scritinv
  zstruct.clambda = llambda
  zstruct.ceta = leta

  mwrfits2, zstruct, zfile, /create, /destroy

  ;; file containing group stats
  ;;mwrfits2, groupstruct, groupfile, /create, /destroy
  mwrfits2, lensum, lensumfile, lhdr, /create, /destroy

  step=oldstep

  ptime, systime(1)-time
  print, 'Finished: '+systime()
  return
END 
