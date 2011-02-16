PRO zobjshear_lambda_setuplens, lensum, rmax, FIRSTLAMBDA, LASTLAMBDA, $
                                stripe, clr, $
                                wz, lclambda, lceta, wlens, DL, angmax, $
                                use_lambda=use_lambda, $
                                clusters=clusters, $
                                hubble=hubble, compcut=compcut,hirata=hirata,$
                                no_sphPolyMask=no_sphPolyMask

  IF n_params() LT 6 THEN BEGIN 
      print,'-Syntax: zobjshear_lambda_setuplens, lensum, rmax, FIRSTLAMBDA, LASTLAMBDA, $'
      print,'                   stripe, clr, $'
      print,'                   wz, lclambda, lceta, wlens, DL, angmax, $'
      print,'                   use_lambda=use_lambda, $'
      print,'                   /clusters, $'
      print,'                   h=h, compcut=compcut,hirata=hirata, $'
      print,'                   no_sphPolyMask=no_sphPolyMask'
      return
  ENDIF 

  IF keyword_set(use_lambda) THEN omegamat=0.3 ELSE omegamat=1.0
  IF n_elements(compcut) EQ 0 THEN compcut = 0.0
  IF n_elements(hubble) EQ 0 THEN hubble=1.0

  max_allow_angle = 6.0         ;degrees

  print,'Max allowed angle: ',max_allow_angle

  nlens = n_elements(lensum)

  ;; get z tag
  wz=getztag(lensum[0])

  ;; get lens clambda,ceta
  getstructlambda, lensum, lclambda, lceta

  ;; Sort by lambda. 
  ;; This is not necessary except for checking systematics along
  ;; stripe...
  sort_stripelambda, lensum, lclambda, lceta

  ;; subscripts for lenses. This will change as we throw
  ;; out lenses
  wlens = lindgen(nlens)

  print,'Using h = ',hubble

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; find DL and thus the maximum search angle
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  angmax = replicate(10000d, nlens)
  DL = replicate(1.e-10, nlens) ; small DL are thrown out in edgecuts
  wtmp = where(lensum.(wz) GT 0)

  DL[wtmp] = angdist_lambda( lensum[wtmp].(wz), $
                             h=hubble, omegamat=omegamat)*1000. 
  angmax = rmax/DL*180.0/!pi     ; angles in degrees

  good = where(angmax LE max_allow_angle)   ; degrees

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; first apply the masks and edge intersections
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; set this if the mask has already been applied
  IF NOT keyword_set(no_sphPolyMask) THEN BEGIN 
      
      IF NOT keyword_set(clusters) THEN BEGIN 
          print,'sphPoly Applying masks'
          print,'and completeness cut: ',compcut
          print
          ;; First apply spectroscopic mask
          apply_sphpoly_mask, lensum[good].ra, lensum[good].dec, bad, good2, $
            compcut=compcut, $
            completeness=completeness
          lensum[good].completeness = completeness
          good = good[good2]
      ENDIF ELSE BEGIN 
          print,'Applying MaxBCG masks'
          print
          MaxBCG_apply_mask, lclambda[good], lceta[good], bad, good2
          good = good[good2]

      ENDELSE 

      nlens1 = n_elements(good)
      print,'---------------------------------------------------------'
      print,'Threw out '+ntostr(nlens-nlens1)+' in lens mask'
      print,'---------------------------------------------------------'
      nlens = nlens1
  ENDIF

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Must apply primary bound cut outside the mask since the mask is
  ;; contains the entire area.  Note, this is only important for stripes
  ;; where spectroscopy goes outside the bounds.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  primary_bound_multi, stripe, bound_array, overall_bound

  good2 = where(lceta[good] LT overall_bound.etamax AND $
                lceta[good] GE overall_bound.etamin, ngood,comp=bad,ncomp=nbad)
  IF ngood EQ 0 THEN BEGIN 
      print,'***************************************'
      message,'No objects passed primary bound cuts'
  ENDIF 

  IF nbad NE 0 THEN lensum[good[bad]].pixelmaskflags = $
    !FLAGS_MASKED + $
    !FLAGS_QUAD1_MASKED+!FLAGS_QUAD2_MASKED+$
    !FLAGS_QUAD3_MASKED+!FLAGS_QUAD4_MASKED

  good = good[good2]
  nlens2 = n_elements(good)
  print,'---------------------------------------------------------'
  print,'Threw out '+ntostr(nlens-nlens2)+' in primary bound cut'
  print,'---------------------------------------------------------'

  nlens = nlens2

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; now apply the pixel mask.  This is now all about how the lenses
  ;; are related to the distribution of sources
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  apply_pixel_mask, lclambda[good], lceta[good], bad, good2, maskflags, $
                    maxangle=angmax[good], /twoquad, /nomissf
  ;;flush,-1,-2
  IF good2[0] EQ -1 THEN BEGIN 
      print,  '*************************************'
      message,'No objects passed masks and edgecuts'
  ENDIF 
  lensum[good].pixelmaskflags = maskflags
  good = good[good2]


  wlens = wlens[good]
  nlens3 = n_elements(wlens)
  print,'----------------------------------------------------------------'
  print,'Threw out '+ntostr(nlens-nlens3)+' in pixel masks and edge cut'
  print,'----------------------------------------------------------------'

  nlens = nlens3

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; mean sigma_crit and DL
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  print,'Calculating mean 1/sigma_crit'
  print
  lensum[wlens].scritinv = $
    sdss_sigma_crit(stripe,clr,lensum[wlens].(wz), $
                    wgood=good2, use_lambda=use_lambda,$
                    hirata=hirata)

  IF good2[0] EQ -1 THEN BEGIN 
      print,  '**************************************'
      message,'No lenses passed sdss_sigma_crit cuts'
  ENDIF 
  wlens = wlens[good2]
  nlens4 = n_elements(wlens)

  print,'---------------------------------------------------------'
  print,'Threw out '+ntostr(nlens-nlens4)+' lenses in sdss_sigma_crit'
  print,'---------------------------------------------------------'
  nlens = nlens4


END 
