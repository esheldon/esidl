PRO conv2d, func, pqlim, ab_limits, xvals, Npts, outfunc

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;  NAME:
;    CONV2D
;
;  PURPOSE:
;    Convolve a symmetric 2d user-defined function 
;    with a 2d user-defined kernel.  These are both evaluated inside the 
;    function 'func' and must access system variables !conv_x0 and 
;    !conv_y0.  pqlim defines the y value to integrate to.  e.g. for a 
;    circle it will integrate to +/-sqrt(R^2 - x^2) 
;    The points of evaluation are the inputs xvals
; 
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: conv2d, func, ab_limits, xvals, outfunc'
      return
  ENDIF 
  
  nx = n_elements(xvals)

; due to symmetry, we can just go along the x-axis
  outfunc = dblarr(nx)

  !conv_y0 = 0.

  FOR ix=0L, nx-1 DO BEGIN 
      !conv_x0 = xvals[ix]

      outfunc[ix]=qgauss2d(func,ab_limits,pqlim, Npts, /double)

      IF ix EQ nx-11 THEN BEGIN
          print,'_10more_',format='(a,$)'
      ENDIF ELSE BEGIN
          print,'.',format='(a,$)'
      ENDELSE 
  ENDFOR 
  print
  return
END
