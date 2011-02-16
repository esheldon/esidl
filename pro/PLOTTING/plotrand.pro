PRO plotrand, x, y, fracuse=fracuse, overplot=overplot, indices=indices, color=color_in, _extra=_extra

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: plotrand, x, y, fracuse=fracuse, overplot=overplot, _extra=_extra'
      print,'Default fracuse = 0.1'
      return
  ENDIF 

  nx=n_elements(x)
  IF n_elements(fracuse) EQ 0 THEN fracuse=1./10.

  ;; plot 1/10 of them, random set
  nuse = long(nx*fracuse)
  
  indices = long(arrscl(randomu(seed,nuse), 0, nx-1, arrmin=0.0, arrmax=1.0))

  newx = x[indices]
  newy = y[indices]

  if n_elements(color_in) ne 0 then color=c2i(color_in)
  IF keyword_set(overplot) THEN BEGIN 
      oplot, newx, newy, color=color, _extra=_extra
  ENDIF ELSE BEGIN 
      plot, newx, newy, color=color, _extra=_extra
  ENDELSE 

END 
