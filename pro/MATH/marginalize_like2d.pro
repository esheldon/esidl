PRO marginalize_like2d, xvals, yvals, likelihood, xlike, ylike

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: marginalize_like2d, xvals, yvals, likelihood, xlike, ylike'
      return
  ENDIF 


  nx = n_elements(xvals)
  ny = n_elements(yvals)

  xlike = dblarr(nx)
  ylike = dblarr(ny)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Integrate over x for each y value
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  npts = 2L*nx
  
  x1 = double(xvals[0])
  x2 = double(xvals[nx-1])
  
  gauleg, x1, x2, npts, XXi, WWi
  
  FOR iy=0L, ny-1 DO BEGIN 
      w=where(likelihood[*,iy] EQ likelihood[*,iy], nw)
      IF nw NE 0 THEN BEGIN 
          funcvals = interpol( likelihood[w, iy], xvals[w], XXi )
          ylike[iy] = total( funcvals*WWi, /double)
      ENDIF 
  ENDFOR 

  ;; normalize
  w = where(ylike NE 0, nw)
  IF nw NE 0 THEN BEGIN 
      ylike = ylike/qgauss(ylike[w], yvals[w], ny*2L)
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Integrate over y for each x value
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  npts = 2L*ny
  
  x1 = double(yvals[0])
  x2 = double(yvals[ny-1])
  
  gauleg, x1, x2, npts, XXi, WWi
  
  FOR ix=0L, nx-1 DO BEGIN 
      w=where(likelihood[ix,*] EQ likelihood[ix,*], nw)
      IF nw NE 0 THEN BEGIN 
          funcvals = interpol( likelihood[ix, w], yvals[w], XXi )
          xlike[ix] = total( funcvals*WWi )
      ENDIF 
  ENDFOR 

  ;; normalize
  w = where(xlike NE 0, nw)
  IF nw NE 0 THEN BEGIN 
      xlike = xlike/qgauss(xlike[w], xvals[w], nx*2L)
  ENDIF 
END 
