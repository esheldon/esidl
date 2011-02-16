PRO oplot_manycolor, x, y, colors, _extra=extra

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: oplot_manycolor, x, y, colors, _extra=extra'
      return
  ENDIF 

  nx=n_elements(x)
  ny=n_elements(y)
  nclrs=n_elements(colors)

  IF (nx NE ny) OR (nx NE nclrs) THEN message,'#x must equal #y and #colors'

  FOR i=0L, nx-1 DO BEGIN 
      oplot,[x[i]], [y[i]], color=colors[i], _extra=extra
  ENDFOR 

END 
