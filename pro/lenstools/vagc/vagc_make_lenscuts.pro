PRO vagc_make_lenscuts, lensum, $
                        rmax, stripes, max_allow_angle, $
                        omegamat, hubble, maxz, compcut, $
                        angMax, DL, wlens, rlrgMask=rlrgMask, $
                        maskFile=maskFile, soFile=soFile

  ;; This is an internal program to be called by, say, vagc_setuplens.pro

  IF n_params() LT 10 THEN BEGIN 

  ENDIF 

  IF keyword_set(rlrgMask) THEN BEGIN 
      maskFile = $
        sdssidl_config('SHAPECORR_DIR')+$
        'masks/pixel_stripe_princeton.mask_simple'
      soFile = '/net/cheops2/home/esheldon/ccode/pixel_masks/apply_maskIDL.so'
  ENDIF

  nlens_init = n_elements(lensum)

  ;; get z tag
  wz=getztag(lensum[0])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find DL and thus the maximum search angle
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  angmax = replicate(10000d, nlens_init)
  DL = replicate(1.e-10, nlens_init) ; small DL are thrown out in edgecuts
  wtmp = where(lensum.(wz) GT 0)

  DL[wtmp] = angdist_lambda( lensum[wtmp].(wz), $
                             h=hubble, omegamat=omegamat)*1000. 
  angmax = rmax/DL*180.0/!pi     ; angles in degrees

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Redshift cuts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wlens = where(angmax LE max_allow_angle AND $
                lensum.(wz) LT maxz, nsubz) ; degrees

  IF nsubz EQ 0 THEN BEGIN 
      print,'***************************************'
      message,'No objects passed redshift cuts'
  ENDIF 

  print,'---------------------------------------------------------'
  print,'Threw out '+ntostr(nlens_init-nsubz)+' in redshift cuts'
  print,'---------------------------------------------------------'  
  nlens = nsubz

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; first apply the completeness cut
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  cgood = where(lensum[wlens].completeness GT compcut, ncomp)

  IF ncomp EQ 0 THEN BEGIN 
      print,'***************************************'
      message,'No objects passed completeness cuts'
  ENDIF 

  print,'---------------------------------------------------------'
  print,'Threw out '+ntostr(nlens-ncomp)+' in completeness cuts'
  print,'---------------------------------------------------------'
  nlens = ncomp
  wlens = wlens[cgood]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;
  ;; Pixel masks and primary bound cuts.
  ;;
  ;; Must apply primary bound cut outside the mask since the mask is
  ;; contains the entire area.  Note, this is only important for stripes
  ;; where spectroscopy goes outside the bounds.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  primary_bound_multi, stripes, bound_array, overall_bound

  pgood = where(lensum[wlens].ceta LT overall_bound.etamax AND $
                lensum[wlens].ceta GE overall_bound.etamin, nprimary, $
                comp=bad,ncomp=nbad)

  IF nprimary EQ 0 THEN BEGIN 
      print,'***************************************'
      message,'No objects passed primary bound cuts'
  ENDIF 

  IF nbad NE 0 THEN lensum[wlens[bad]].pixelmaskflags = $
    !FLAGS_MASKED + $
    !FLAGS_QUAD1_MASKED+!FLAGS_QUAD2_MASKED+$
    !FLAGS_QUAD3_MASKED+!FLAGS_QUAD4_MASKED


  print,'---------------------------------------------------------'
  print,'Threw out '+ntostr(nlens-nprimary)+' in primary bound cut'
  print,'---------------------------------------------------------'
  nlens = nprimary
  wlens = wlens[pgood]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; now apply the pixel mask.  This is now all about how the lenses
  ;; are related to the distribution of sources
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  apply_pixel_mask, $
    lensum[wlens].clambda, lensum[wlens].ceta, masked, unmasked, maskFlags, $
    maxangle=angmax[wlens], /twoquad, /nomissf, $
    soFile=soFile, maskFile=maskFile

  ;;flush,-1,-2
  IF unmasked[0] EQ -1 THEN BEGIN 
      print,  '*************************************'
      message,'No objects passed masks and edgecuts'
  ENDIF 
  lensum[wlens].pixelMaskFlags = maskFlags

  nunmasked = n_elements(unmasked)

  print,'----------------------------------------------------------------'
  print,'Threw out '+ntostr(nlens-nunmasked)+' in pixel masks and edge cut'
  print,'----------------------------------------------------------------'

  wlens = wlens[unmasked]
  nlens = nunmasked

  ;; No need for this really. If so, will have to re-add clr, /hirata stuff?
return
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; mean sigma_crit and DL
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  print
  print,'Calculating mean 1/sigma_crit'
  print
  lensum[wlens].scritinv = $
    sdss_sigma_crit(stripes,clr,lensum[wlens].(wz), $
                    wgood=scritgood, use_lambda=use_lambda,$
                    hirata=hirata)

  IF scritgood[0] EQ -1 THEN BEGIN 
      print,  '**************************************'
      message,'No lenses passed sdss_sigma_crit cuts'
  ENDIF 

  nscritgood = n_elements(scritgood)

  print,'---------------------------------------------------------'
  print,'Threw out '+ntostr(nlens-nscritgood)+' lenses in sdss_sigma_crit'
  print,'---------------------------------------------------------'
 
  wlens = wlens[scritgood]

END 
