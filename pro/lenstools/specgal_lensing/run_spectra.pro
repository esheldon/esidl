PRO run_spectra, stripes, clr, nbin_OR_binsize, $
                 indir=indir, $
                 outdir=outdir, $
                 $
                 rmin=rmin, $
                 rmax=rmax, $
                 use_lambda=use_lambda, $
                 hirata=hirata, $
                 photoz=photoz, $
                 recorr=recorr,$
                 $
                 logbin=logbin, $
                 $
                 depth=depth, $
                 $
                 extno=extno, $
                 $
                 lrg_lenses=lrg_lenses, $
                 lrg_sources=lrg_sources, $
                 rlrg_sources=rlrg_sources, $
                 newPhotoZ=newPhotoZ, $
                 $
                 pixelMaskFlags=pixelMaskFlags, $
                 rlrgMask=rlrgMask, $
                 $
                 print2ascii=print2ascii, $
                 $
                 scat=scat, revind=revind, $
                 lenscat=lenscat, $
                 $
                 randnum=randnum, $
                 noTestQuad=noTestQuad, $
                 $
                 callCPP=callCPP

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax:  run_spectra, stripes, clr, nbin_OR_binsize, $'
      print,'    indir=indir, $'
      print,'    outdir=outdir, $'
      print,'    $'
      print,'    rmin=rmin, $'
      print,'    rmax=rmax, $'
      print,'    /use_lambda, $'
      print,'    /hirata, $'
      print,'    /photoz, $'
      print,'    /recorr,$'
      print,'    $'
      print,'    /logbin, $'
      print,'    $'
      print,'    depth=depth, $'
      print,'    $'
      print,'    extno=extno, $'
      print,'    $'
      print,'    /lrg_lenses, $'
      print,'    /lrg_sources, $'
      print,'    /rlrg_sources, $'
      print,'    /newPhotoZ, $'
      print,'    $'
      print,'    /pixelMaskFlags, $'
      print,'    /rlrgMask, $'
      print,'    $'
      print,'    /print2ascii, $'
      print,'    $'
      print,'    scat=scat, revind=revind, $'
      print,'    lenscat=lenscat, $'
      print,'    $'
      print,'    randnum=randnum, $'
      print,'    /noTestQuad, $'
      print,'    /callCPP'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  colors=['u','g','r','i','z']

  sdssidl_setup,/silent
  setup_mystuff

  IF n_elements(rmin)       EQ 0 THEN rmin       = 20.
  IF n_elements(rmax)       EQ 0 THEN rmax       = 1020.
  IF n_elements(use_lambda) EQ 0 THEN use_lambda = 1
  IF n_elements(hirata)     EQ 0 THEN hirata     = 1
  IF n_elements(photoz)     EQ 0 THEN photoz     = 1
  IF n_elements(recorr)     EQ 0 THEN recorr     = 1
  
  IF n_elements(outdir) EQ 0 THEN $
    outdir = esheldon_config("lensout_dir")+'stripe'+$
    stripearr2string(stripes)+'/'

  IF NOT fexist(outdir) THEN BEGIN
      spawn,'mkdir '+outdir
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; lens catalog
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(lenscat) EQ 0 THEN BEGIN 
      vagc_lensinput_name, stripes, rmin, rmax, name, $
        randnum=randnum, rlrgMask=rlrgMask
      IF NOT fexist(name) THEN BEGIN 
          message,'Lens input file does not exist: '+name
      ENDIF 
      
      ;; order is important since we are sending this to the C++ code
      columns = ['ra','dec','clambda','ceta',$
                 'z',$
                 'zindex','index',$
                 'scritinv','DL','angmax',$
                 'pixelmaskflags',$
                 'totpairs',$
                 'ie',$
                 'weight',$
                 $
                 'angsum',$
                 'rsum',$
                 'rmin_act',$
                 'rmax_act',$
                 $
                 'sigma',$
                 'sigmaerr',$
                 'sigerrsum',$
                 'orthosig',$
                 'orthosigerr', $
                 'orthosigerrsum',$
                 $
                 'sshsum',$
                 'wsum',$
                 'owsum',$
                 'wsum_ssh',$
                 'npair']

      print
      print,'Reading lens input file: ',name
      lenscat = mrdfits(name, 1, columns=columns)

      IF n_elements(tag_names(lenscat[0])) NE n_elements(columns) THEN BEGIN 
          message,'Error reading struct columns'
      ENDIF 
  ENDIF 

  lenscat.scritinv = $
    sdss_sigma_crit(stripes,clr,lenscat.z, $
                    wgood=scritgood, use_lambda=use_lambda,$
                    hirata=hirata, rlrg_sources=rlrg_sources)
  lenscat = lenscat[scritgood]

  ngal = n_elements(lenscat)
  print
  print,'Lenses: ',ngal

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; get source catalog. 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(rlrg_sources) THEN BEGIN 
      sfile = '~/Rachel/LRGscat.fit'
      
      print
      print,'Reading file: ',sfile
      scat = mrdfits(sfile,1)

      IF keyword_set(newPhotoZ) THEN BEGIN 
          scat.photoz_z = scat.photoz_z_lrg
          deltaFuncPhotoZ = 1
      ENDIF 
  ENDIF ELSE IF keyword_set(callCPP) THEN BEGIN 

      IF n_elements(scat) EQ 0 OR n_elements(revind) EQ 0 THEN BEGIN 

          file = $
            sdssidl_config('shapecorr_dir') + 'combined/'+$
            'stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_76_82_86_srcgal_gri_h_nzCuts_specGal.fit'
          
          print
          print,'Reading scat: ',file
          scat = mrdfits(file,1,hdr)
          print,'Reading revind'
          revind = mrdfits(file,2)

          htmDepth = sxpar(hdr, 'HTMDEPTH')

          IF n_elements(depth) NE 0 THEN BEGIN 
              IF depth NE htmDepth THEN BEGIN 
                  
                  csurvey2eq, scat.clambda, scat.ceta, ra, dec
                  print
                  print,'Looking up leafids for depth = ',depth
                  htmLookupRadec, ra, dec, depth, leafids
                  scat.leafid = leafids

                  print,'Sorting scat'
                  s = sort(scat.leafid)
                  scat = scat[s]

                  print,'Histograming leafids'
                  minid = min(scat.leafid, max=maxid)
                  leaf_hist = histogram(scat.leafid, min=minid, max=maxid, $
                                        reverse_indices=revind)

              ENDIF 
          ENDIF 
          

      ENDIF 

  ENDIF ELSE IF n_elements(scat) EQ 0 THEN BEGIN 
      IF n_elements(scat) EQ 0 OR n_elements(revind) EQ 0 THEN BEGIN 

          file = $
            sdssidl_config('shapecorr_dir') + 'combined/'+$
            'stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_76_82_86_srcgal_gri_h_nzCuts_specGal.fit'
          
          print
          print,'Reading scat: ',file
          scat = mrdfits(file,1)
          print,'Reading revind'
          revind = mrdfits(file,2)
      ENDIF 

;      columns = ['clambda','ceta','e1_recorr','e2_recorr',$
;                 'e1e1err','e1e2err','e2e2err','photoz_z','photoz_zerr']

      ;; Not doing pixel mask cuts here. 
;      get_scat, stripes, clr, scat, /hirata, $
;        pixelMaskFlags=pixelMaskFlags, rlrgMask=rlrgMask, /nzcuts, $
;        columns=columns
  ENDIF 

;  IF keyword_set(lrg_lenses) THEN BEGIN 
;      addstr='lrg_zgal_gal'
;      raddstr = 'lrg'
;  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; measure shear around lenses
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(randnum) NE 0 THEN BEGIN 
      addstr = 'zrand'+ntostr(randnum)
  ENDIF 

  step = long(300*(ngal/30000.))
  zobjshear_lambda, $
    stripes, lenscat, scat, clr, rmin, rmax, nbin_or_binsize, $
    callCPP=callCPP, revind=revind, $
    step=step, use_lambda=use_lambda, $  
    sumfile=sumfile, lensumfile=lensumfile, $
    outdir=outdir,$
    logbin=logbin,$
    depth=depth, photoz=photoz, recorr=recorr,$
    extno=extno, $
    hirata=hirata, compcut=compcut, $
    print2ascii=print2ascii, addstr=addstr, $
    lrg_sources=lrg_sources, rlrg_sources=rlrg_sources, $
    noTestQuad=noTestQuad, deltaFuncPhotoZ=deltaFuncPhotoZ

return

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; now around random points
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  delvarx,scat
  get_scat, stripes, clr, scat, indir=scatdir, hirata=hirata, columns=columns

jump:

  delvarx,lenscat

  tmp=mrdfits(zfile,1)
  zrand = tmp.z

  numrand = ndup*n_elements(zrand)
  step = long(300*(numrand/30000.))

  zrandshear_lambda, $
    stripes, scat, clr, rmin, rmax, nbin_or_binsize, zrand, ndup, $
    step=step,$
    addstr=raddstr,$
    outdir=outdir, use_lambda=use_lambda,zfile=rzfile,$
    logbin=logbin,$
    depth=depth, photoz=photoz, recorr=recorr,$
    extno=extno, $
    hirata=hirata, compcut=compcut, $
    print2ascii=print2ascii, lrg_sources=lrg_sources

  delvarx, scat, tmp, zrand
  
  return
END 
