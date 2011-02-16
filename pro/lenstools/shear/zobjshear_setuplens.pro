PRO zobjshear_setuplens, lensum, rmax, FIRSTRA, LASTRA, $
                         stripe, clr, $
                         wz, lra, ldec, wlens, DL, angmax, $
                         tsgals=tsgals, maskdir=maskdir, $
                         use_lambda=use_lambda, commonsrc=commonsrc, $
                         hubble=hubble,$
                         issouth=issouth     

  IF n_params() LT 6 THEN BEGIN 
      print,'-Syntax: zobjshear_setuplens, lensum, rmax, FIRSTRA, LASTRA, stripe, clr, $'
      print,'     wz, lra, ldec, wlens, DL, angmax, $'
      print,'     tsgals=tsgals, maskdir=maskdir, $'
      print,'     use_lambda=use_lambda, commonsrc=commonsrc, hubble=hubble'
      return
  ENDIF 

  ninit = n_elements(lensum)

  ;; get z tag
  wz=getztag(lensum[0])
  ;; get lens ra,dec
  getstructlambda, lensum, lra, ldec, /radec

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Sort by ra. If it is a southern stripe, deal with fact that it
  ;; crosses [0,360] mark (possibly)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sort_stripera, lensum, lra, ldec, issouth

  ;; subscripts for lenses. This will change as we throw
  ;; out lenses
  wlens = lindgen(ninit)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; sigma_crit and DL
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  lensum.scritinv = sdss_sigma_crit(stripe, clr, lensum.(wz[0]), $
                                    wgood=wgood, use_lambda=use_lambda, $
                                    commonsrc=commonsrc)

  print,'Using h = ',h
  print,'Calculating 1/sigma_crit'
  print

  ;;Convert from Mpc to kpc
  IF NOT keyword_set(use_lambda) THEN omegamat=1.0 ELSE omegamat=0.3
  DL = angdist_lambda( lensum.(wz[0]), h=hubble, omegamat=omegamat)*1000. 
  angmax = rmax/DL*180./!pi     ;angles in degrees

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; don't want angle to be larger than 1/rfac of a stripe (until we use 
  ;; multiple stripes)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  rfac = 1.0
  max_allowed_angle = 2.5/rfac
  wgood2 = where(angmax[wgood] LE max_allowed_angle)
  wgood = wgood[wgood2]
  wlens = wlens[wgood]

  nlens1 = n_elements(wlens)
  print,'Threw out ',ntostr(ninit - nlens1),' lenses because too deep or close'

  ;; apply mask (rotated frame)

  read_stripe_mask, stripe, mask, tsgals=tsgals, indir=maskdir
  IF issouth THEN BEGIN
      rotate_ra, lra
      rotate_ra, FIRSTRA
      rotate_ra, LASTRA
  ENDIF 

  ;; apply ra cut
  wtmp = where( ((lra[wlens] - angmax[wlens]) GE FIRSTRA ) AND $
                ((lra[wlens] + angmax[wlens]) LE LASTRA ), nlens)
  IF nlens EQ 0 THEN message,'No objects passed ra cut'
  wlens = wlens[wtmp]
  print,'Threw out ',ntostr(nlens1-nlens),' lenses in ra cut'

  IF issouth THEN BEGIN
      rotate_ra, lra
      rotate_ra, FIRSTRA
      rotate_ra, LASTRA
  ENDIF 


END 
