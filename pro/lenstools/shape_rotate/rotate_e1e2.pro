PRO rotate_e1e2, angle, e1_in, e2_in, e1_out, e2_out

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: rotate_e1e2, angle, e1_in, e2_in [, e1_out, e2_out]'
      print,'angle in radians'
      return
  ENDIF 

  ;; angle should either be same size as e1,e2 or a single number
  nang=n_elements(angle)
  ne1=n_elements(e1_in)
  ne2=n_elements(e2_in)

  IF ne1 NE ne2 THEN BEGIN 
      message,'e1_in and e2_in must be same size',/inf
      return
  ENDIF 

  bad=1
  IF nang EQ 1 THEN bad=0
  IF bad THEN BEGIN
      IF (nang NE ne1) THEN BEGIN 
          message,'If angle is an array, it must be same size as e1_in,e2_in',/inf
          return
      ENDIF ELSE bad=0 
  ENDIF 

  cos2a = cos(2.*angle)
  sin2a = sin(2.*angle)
  
  e1_out =  e1_in*cos2a + e2_in*sin2a
  e2_out = -e1_in*sin2a + e2_in*cos2a

  return
END 
  
