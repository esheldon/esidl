FUNCTION shearsis, sigma, zsource, zlens, angle, kpc=kpc

  IF n_params() LT 3 THEN BEGIN
      print,'Syntax:  s = shearsis(sigma, zsource, zlens, angle, kpc=kpc)'
      print,'Sigma in km/s.   Angle in arcsec.   Assumes matter only universe'
      print,'kpc overrides angle.  Forces physical units.'
      return,-1
  ENDIF 

  IF n_elements(kpc) EQ 0 THEN BEGIN 

      ratio = angdist(zsource,zlens,h=1.)/angdist(zsource,h=1.)
      return, 1.4*(sigma/220.)^2*ratio/angle/2.0

  ENDIF ELSE BEGIN
      
      c=3.e5                    ;km/s
      
      D = angdist(zlens,h=1.)*angdist(zsource,zlens,h=1.)/angdist(zsource,h=1.)
                                ; Mpc
      R = kpc/1000.             ; Mpc
      return, 2.*!pi*(sigma/c)^2*D/R
  ENDELSE 

END 
