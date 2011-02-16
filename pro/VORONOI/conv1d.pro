PRO conv1d, func, Alim, Blim, xvals, Npts, outfunc

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;  NAME:
;    CONV1D
;
;  PURPOSE:
; 
; 
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF n_params() LT 4 THEN BEGIN
      print,'-Syntax: '
      return
  ENDIF 

  nx = n_elements(xvals)
  outfunc = dblarr(nx)

  FOR ix=0L, nx-1 DO BEGIN 
      print,'.',format='(a,$)'
      !conv_x0 = xvals[ix]

      outfunc[ix]=QGAUSS(func, Alim, Blim, Npts)

  ENDFOR 
  
  return
END
