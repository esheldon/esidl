FUNCTION sigmasis_trunc, sigma, cutoff, radius, core=core, contrast=contrast

  IF n_params() LT 3 THEN BEGIN
      print,'Syntax:  s = sigmasis_trunc(sigma, cutoff, zsource, zlens, radius, core=core, contrast=contrast)'
      print,'if /core, then add 1kpc core'
      print,'Sigma in km/s.   Radius/cutoff in kpc.'
      return,-1
  ENDIF 

  IF keyword_set(core) THEN docore = 1 ELSE docore = 0

  IF NOT keyword_set(contrast) THEN BEGIN 
      IF docore THEN BEGIN
          term =  1./sqrt(1. + radius^2)
      ENDIF ELSE BEGIN
          term =  1./radius
      ENDELSE 
      return, 3352.40*(sigma/170.)^2*( term - 1./sqrt(radius^2 + cutoff^2 ) ) ;Msolar/pc^2
  ENDIF ELSE BEGIN 
      IF docore THEN BEGIN
          term1 =  1./sqrt(1. + radius^2)
          term2 =  2.*cutoff*term1^2
          term3 = -1./sqrt( radius^2 + cutoff^2 )
          term4 =  2.*cutoff^2*term1^2*term3
      ENDIF ELSE BEGIN
          term1 =  1./radius
          term2 =  2.*cutoff*term1^2
          term3 = -1./sqrt( radius^2 + cutoff^2 )
          term4 =  2.*cutoff^2*term1^2*term3
      ENDELSE 
      return, 3352.40*(sigma/170.)^2*( term1+term2+term3+term4 ) ;Msolar/pc^2
  ENDELSE 


END 
      
