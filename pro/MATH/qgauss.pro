pro _qgauss_checkc

  common gauss1dblk, use_cproc, xxi, wwi
  if n_elements(use_c_write) eq 0 then begin 
      proNames  = routine_info(/system)
      w = where(proNames EQ 'GAULEG',nw)
      if nw eq 0 then use_cproc = 0 else use_cproc = 1
  endif 

end

PRO _qgauss_gauleg, npts, silent=silent

  common gauss1dblk, use_cproc, xxi, wwi

  if n_elements(use_cproc) eq 0 then _qgauss_checkc

  ;; only need to do this if npts changes
  IF npts NE n_elements(WWi) THEN BEGIN 
      IF NOT keyword_set(silent) THEN BEGIN 
          print,'QGAUSS: Calculating weights and abscissa for npoints =  ',npts
      ENDIF 
      if use_cproc then begin
          gauleg, -1., 1., npts, XXi, WWi
      endif else begin
          gauleg_idl, -1., 1., npts, XXi, WWi
      endelse
  ENDIF 

END 

FUNCTION _qgauss_func, func, x1, x2, npts, silent=silent

  common gauss1dblk, use_cproc, xxi, wwi

  _qgauss_gauleg, npts, silent=silent

  ;; these needed for coordinate transformation
  f1 = (x2-x1)/2.0
  f2 = (x2+x1)/2.0
  
  ;; interpolate
  xvals = XXi*f1 + f2
  return, total(CALL_FUNCTION(func, xvals)*WWi)*f1

END 

FUNCTION _qgauss_data, y, x, npts, silent=silent

  common gauss1dblk, use_cproc, xxi, wwi

  _qgauss_gauleg, npts, silent=silent
  
  ;; Note x must be ordered
  x1 = x[0]
  x2 = x[n_elements(x)-1]
  
  ;; these needed for coordinate transformation
  f1 = (x2-x1)/2.0
  f2 = (x2+x1)/2.0
  
  ;; interpolate
  xvals = XXi*f1 + f2
  yvals = interpol(y, x, xvals)
  return, total(yvals*WWi)*f1


END 

FUNCTION qgauss, func, x1, x2, npts, silent=silent, $
                 deja_vu=deja_vu ; ignored

  on_error, 2
  np = n_params()
  IF np NE 3 AND np NE 4 THEN BEGIN
      print,'-Syntax: '
      print,'         result=qgauss(y, x, npts, /silent)'
      print,' OR'
      print,'         result=qgauss(funcname, x1, x2, npts, /silent)'
      print,'            Where funcname is a string name of a function '
      print,'            taking only abscissa'
      print
      message,'Halting'
  ENDIF 
  
  IF size(func, /type) EQ 7 AND np EQ 4 THEN BEGIN 
      return, _qgauss_func(func, x1, x2, npts, silent=silent)
  ENDIF ELSE BEGIN 
      return, _qgauss_data(func, x1, x2, silent=silent)
  ENDELSE 

END 


