PRO rotate_e1e2error, angle, e1, e2, e1err, e1e2err, e2err, $
                      e1_out, e2_out, e1err_out, e1e2err_out, e2err_out

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: rotate_e1e2error, angle, e1_in, e2_in, e1err_in, e1e2err_in, e2err_in, $'
      print,'                                  e1_out, e2_out, e1err_out, e1e2err_out, e2err_out'
      print,'angle in radians'
      return
  ENDIF 

  ;; angle should either be same size as e1,e2 or a single number
  nang=n_elements(angle)
  ne1=n_elements(e1)
  ne2=n_elements(e2)
  ne1err = n_elements(e1err)
  ne2err = n_elements(e2err)
  ne1e2err = n_elements(e1e2err)

  IF ((ne1 NE ne2) OR $
      (ne1 NE ne1err) OR $
      (ne1 NE ne2err) OR $
      (ne1 NE ne1e2err)) THEN BEGIN 
      message,'e* arrays must all be same size'
  ENDIF 

  diff=1
  IF nang EQ 1 THEN diff=0
  IF diff THEN BEGIN
      IF (nang NE ne1) THEN BEGIN 
          message,'If angle is an array, it must be same size as e1,e2'
      ENDIF
  ENDIF 

  cos2a = cos(2.*angle)
  sin2a = sin(2.*angle)
  
  e1_out = fltarr(ne1)
  e2_out = e1_out
  e1err_out = e1_out
  e2err_out = e2_out
  e1e2err_out = e1_out

  IF nang EQ 1 THEN BEGIN
      R =   [ [ cos2a[0], sin2a[0]], $
              [-sin2a[0], cos2a[0]] ]
  ENDIF 

  FOR i=0L, ne1-1 DO BEGIN 

      e = [e1[i], e2[i]]
      
      e1e2err2 = e1e2err[i]^2*sign(e1e2err[i])
      var = [ [e1err[i]^2, e1e2err2  ],     $
              [e1e2err2,   e2err[i]^2] ]
      
      IF nang GT 1 THEN BEGIN 
          R =   [ [ cos2a[i], sin2a[i]], $
                  [-sin2a[i], cos2a[i]] ]
      ENDIF 

      eprime = R##e
      varprime = R##(var##transpose(R))
            
      e1_out[i] = eprime[0]
      e2_out[i] = eprime[1]
      e1err_out[i] = sqrt(varprime[0,0])
      e2err_out[i] = sqrt(varprime[1,1])
      e1e2err_out[i] = sqrt(abs(varprime[0,1]))*sign(varprime[0,1])
      
  ENDFOR 


END 
