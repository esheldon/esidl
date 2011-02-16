
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
;    zobjshear_lambda, stripes, lenscat, scat, clr, rmin, rmax, binsize, $
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
;    stripes: stripes number in integer form.
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
;    /hirata: useing the chris hirata corrections, modify name
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


PRO zobjshear_lambda_callcpp, lensum, scat, revind, scinvStruct, par_struct

  sofile = '/net/cheops2/home/esheldon/ccode/objShear/bin/objShearIDL.so'
  entry = 'main'

  tt = systime(1)
  tmp = call_external(sofile, entry, $
                      $
                      lensum, scat, revind, scinvStruct, par_struct, $
                      value=[0B, 0B, 0B, 0B, 0B], /unload)

  print,'Time for C++ code: '
  ptime,systime(1) - tt

END 

PRO zobjshear_lambda, stripes, lensum, scat, clr, rmin, rmax, $
                      callCPP=callCPP,revind=revind, $
                      nbin_or_binsize, $
                      use_lambda=use_lambda, $
                      step=step, addstr=addstr, $
                      outdir=outdir, $
                      wgood=wgood, $
                      check=check, $
                      datfile=datfile, sumfile=sumfile, zfile=zfile, $
                      lensumfile=lensumfile, $
                      maxe=maxe, logbin=logbin, $
                      comoving=comoving, $
                      depth=depth, recorr=recorr, $
                      extno=extno, compcut=compcut, combined=combined,$
                      hirata=hirata, print2ascii=print2ascii, $
                      clusters=clusters, $
                      lrg_sources=lrg_sources, rlrg_sources=rlrg_sources, $
                      photoz=photoz, deltaFuncPhotoZ=deltaFuncPhotoZ, $
                      noTestQuad=noTestQuad, checkRidgeLine=checkRidgeLine

  IF n_params() LT 7 THEN BEGIN
      print,'-Syntax:'
      return
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  setup_mystuff

  time = systime(1)

  IF n_elements(maxe) EQ 0 THEN maxe = .2
  IF n_elements(depth) EQ 0 THEN depth=10

  ;; Size by which to group things  300 was shown to be fastest in certain
  ;; circumstances
  IF n_elements(step) EQ 0 THEN step = 300L
  oldstep = step
  IF n_elements(addstr) EQ 0 THEN addstr = 'zgal_gal'
  IF n_elements(outdir) EQ 0 THEN $
    outdir = esheldon_config("lensout_dir")+'tmp/'
  IF NOT keyword_set(check) THEN check=0

  stripestr = '_stripe'+stripearr2string(stripes)
  clr_str = clrarr2string(clr)

  IF keyword_set(recorr) THEN recorrstr='_recorr' ELSE recorrstr=''

  ;; hirata is default
  IF n_elements(hirata) EQ 0 THEN hirata=1
  IF keyword_set(hirata) THEN hiratastr = '_h' ELSE hiratastr=''

  ;; use_lambda default
  IF n_elements(use_lambda) EQ 0 THEN use_lambda=1
  IF keyword_set(use_lambda) THEN omega_m=0.3 ELSE omega_m=1.0
 
  IF keyword_set(lrg_sources) THEN BEGIN 
      lrgstr = '_lrg' 
  ENDIF ELSE BEGIN 
      IF keyword_set(rlrg_sources) THEN lrgstr = '_rlrg' ELSE lrgstr=''
  ENDELSE 

  nlens = n_elements(lensum)
  nscat = n_elements(scat)

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
      logbin=0
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
  hubble=1.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up output files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; if extno is sent, then does not check for existence of files
  prename = $
    outdir+addstr+stripestr+'_'+clr_str+lrgstr+recorrstr+hiratastr
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
; if /photoz, get scinvStruct to interpolate
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; this already includes the prior
  IF keyword_set(photoz) AND NOT keyword_set(deltaFuncPhotoZ) THEN BEGIN 
      get_meanscinv, stripes, clr, scinvStruct, use_lambda=use_lambda, $
        hirata=hirata, $
        lrg_sources=lrg_sources, rlrg_sources=rlrg_sources

      interpPhotoz = 1
  ENDIF ELSE interpPhotoz = 0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up parameter struct
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  minLeafId = min(scat.leafid, max=maxLeafId)
  par_struct = create_struct('h',             float(hubble), $
                             'omega_m',       float(omega_m), $
                             $
                             'interpPhotoz',  fix(interpPhotoz),$
                             $
                             'logbin',        fix(logbin),$
                             $
                             'nbin',          fix(nbin),$
                             $
                             'binsize',       float(binsize),$
                             $
                             'rmin',          float(rmin), $
                             'rmax',          float(rmax), $
                             $
                             'comoving',      keyword_set(comoving), $
                             $
                             'nlens',         long(nlens), $
                             'nsource',       long(nsource), $
                             $
                             'depth',         fix(depth), $
                             'minLeafId',     long(minLeafId), $
                             'maxLeafId',     long(maxLeafId) )


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; plot out the redshifts and positions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  begplot,name=psfile

  !p.multi=[0,0,2]

  wz = getztag(lensum[0])
  plothist, lensum.(wz), bin=0.01, xtitle='Z!DLens!N', ytitle='N(Z)'

  plot, lensum.clambda, lensum.ceta, psym=3, /ynozero, $
        xtitle = !csym.lambda+'!Dc!N', ytitle = !csym.eta+'!Dc!N'
  !p.multi=0
  endplot 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Estimate shear around all of these lenses
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(callCPP) THEN BEGIN 

      help,revind
      IF n_elements(revind) EQ 0 THEN message,'you must enter revind'
;;      zobjshear_lambda_callcpp, lensum, scat, revind, scinvStruct, par_struct
;stop
      objshear, lensum, scat, revind, scinvStruct, par_struct
;;stop
      ;; make cuts based on ie and weight
      wlens = where(lensum.weight GT 0)
      mm = 3.0/sqrt(lensum[wlens].totPairs)

      mmm = mm > maxe
      wlens2 = where(lensum[wlens].ie LE mmm, lensused)
      wlens = wlens[wlens2]

  ENDIF ELSE BEGIN 

      htmTime = systime(1)
      zobjshear_htm_looplens, depth, lensum, scat,$
        step, maxe, $
        rmin, rmax, nbin_or_binsize, $
        groupstruct, indices, $
        logbin=logbin, $
        recorr=recorr,$
        scinvStruct=scinvStruct, $
        outlun=outlun, $
        deltaFuncPhotoZ=deltaFuncPhotoZ, $
        noTestQuad=noTestQuad, checkRidgeLine=checkRidgeLine

      IF keyword_set(print2ascii) THEN free_lun, outlun

      print,'Done with looplens'
      ptime, systime(1)-htmTime

      wlens = where(indices NE -1, lensused)

  ENDELSE 

  ;; Remove unused or bad lenses
  
  print
  echo,'Finally used '+ntostr(lensused)+'/'+ntostr(nlens)+' Lenses',$
    color='green',/bold

  ;; delete the scat to free memory
  delvarx, scat
  lensum = lensum[wlens]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Find averages from sums
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



  lensw = lensum.weight
  lenswsum = total( lensw )
  zsum = total(lensum.(wz)*lensw)
  zmean = zsum/lenswsum
  print
  print,'zmean = ',zmean

  combine_zlensum, lensum, binsize, rmin, rmax, hubble, shstruct, $
    compcut=compcut, depth=depth, comoving=comoving, logbin=logbin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Ouput the data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; Shear file for these lenses
  shhdr = zshhdr(shstruct)  
  lhdr  = zsumhdr(shstruct)
  print
  print,'Writing mean file: ',datfile
  mwrfits, shstruct, datfile, shhdr, /create

  ;; File containing ra,dec,redshift for each used lens
  zs = create_struct('z', 0., $
                     'scritinv', 0., $
                     'clambda', double(0.), $
                     'ceta', double(0.) )

  zstruct          = replicate(zs, lensused)
  zstruct.z        = lensum.(wz)
  zstruct.scritinv = lensum.scritinv
  zstruct.clambda  = lensum.clambda
  zstruct.ceta     = lensum.ceta

  print,'Writing z file: ',zfile
  mwrfits, zstruct, zfile, /create

  ;; file containing group stats
  print,'Writing group file: ',datfile
  mwrfits, groupstruct, groupfile, /create

  ;; File with all the used lenses
  print,'Writing lensum file: ',lensumfile
  mwrfits2, lensum, lensumfile, lhdr, /create, /destroy

  step=oldstep
  ptime, systime(1)-time
  print, 'Finished: '+systime()
  return
END 
