FUNCTION binzvol, z, dz, omega=omega, h=h

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; Name:  binzvol
;
; Purpose: calculate the volume in a thin (in z) shell.  Thin is 
;          defined as dz << z
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF n_params() EQ 0 THEN BEGIN
    print,'-Syntax: binzvol, z, dz, v, imega=omega, h=h'
    print,' Returns volume in cubic Mpc'
    return,0.
  ENDIF 

  IF NOT keyword_set(omega) THEN omega = 1.0
  IF NOT keyword_set(h) THEN h = .7

  ;Define some parameters
  c = 3.0e5 ;km/s
  H0 = 100.*h

  ;find angular diameter distance to z
  dang1 = angdist(z)
  dang2 =  angdist(z+dz)
  delta = dang2-dang1

  ; approximate volume element
  v = 4./3.*!dpi*3.*delta*dang1^2

  return,v
END
          
