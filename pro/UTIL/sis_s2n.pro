FUNCTION sis_s2n, rad1, rad2, sigma, zs, zl, Nlens, arcsec=arcsec

  ;; Finds the s2n of shear in a bin R1 and R2 wide in physical coordinates
  ;; or angular

  IF n_params() LT 4 THEN BEGIN
      print,'-Syntax: result=sis_s2n(R1, R2, sigma, zs, zl, Nlens, arcsec=arcsec'
      print,' sigma in km/s'
      return,-1
  ENDIF 

  IF keyword_set(arcsec) THEN BEGIN 
      conv = 1./3600.*!pi/180.
      R1 = rad1*conv              ;rad
      R2 = rad2*conv              ;rad

      R1 = R1*angdist(zl)       ;Mpc
      R2 = R2*angdist(zl)
  ENDIF ELSE BEGIN
      R1 = rad1
      R2 = rad2
  ENDELSE 

  Area = !pi*(R2^2 - R1^2)
  Rm = gmean(R1, R2, 2.)        ;geometric mean
  
  Am = !pi*Rm^2

  
  ratio=angdist(zs,zl)/angdist(zs)
  fac=8.5e-7
  lam=2.                        ;# per square arcminute
  sigshape=.55

  s2n=fac*sigma^2*ratio*sqrt(lam*Area/Am )*sqrt(Nlens)/sigshape

  return,s2n

END 
