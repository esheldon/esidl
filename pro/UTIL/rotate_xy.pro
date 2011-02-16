PRO rotate_xy, angle, x_in, y_in, x_out, y_out

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: rotate_xy, angle, x_in, y_in [, x_out, y_out]'
      print,'angle in radians'
      return
  ENDIF 

  ;; angle should either be same size as x,y or a single number
  nang=n_elements(angle)
  nx=n_elements(x_in)
  ny=n_elements(y_in)

  IF nx NE ny THEN BEGIN 
      message,'x_in and y_in must be same size',/inf
      return
  ENDIF 

  bad=1
  IF nang EQ 1 THEN bad=0
  IF bad THEN BEGIN
      IF (nang NE nx) THEN BEGIN 
          message,'If angle is an array, it must be same size as x_in,y_in',/inf
          return
      ENDIF ELSE bad=0 
  ENDIF 

  cos2a = cos(2.*angle)
  sin2a = sin(2.*angle)
  
  x_out =  x_in*cos2a + y_in*sin2a
  y_out = -x_in*sin2a + y_in*cos2a

  return
END 
  
