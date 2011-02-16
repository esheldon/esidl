function gaussprob, xi, mean, sigma, DOUBLE=double, minx=minx, useind=useind

  Nparms = n_params()
  if Nparms LT 3 then begin
        print,'Syntax - y = GAUSSPROB( xi, mean, sigma, /DOUBLE ])'
        return, -1
  endif

  npts = N_elements(xi) 
;  mingauss = (machar(DOUBLE=double)).xmin
mingauss = 1.e-10
  gauss = replicate(mingauss, npts )
  IF keyword_set(double) THEN BEGIN 
      ;;gauss = dblarr( npts ) 
      mean = double(mean)
      sigma = double(sigma)
  ENDIF ELSE BEGIN 
      ;;gauss = fltarr( npts )
      mean = float(mean)
      sigma = float(sigma)
  ENDELSE 
  
  z = ( xi - mean )/(sigma*sqrt(2.))
  zz = z*z
  
  norm = 1./sqrt(2.*!dpi)/sigma

; Get smallest value expressible on computer.   Set lower values to 0 to avoid
; floating underflow
  minexp = alog(mingauss)     
  
  w = where( zz LT -minexp, nw )
  if (nw GT 0) then gauss[w] = norm*exp( -zz[w] )
  
  IF n_elements(minx) NE 0 THEN BEGIN 
      useind=where(xi GE minx, nuseind, comp=comp, ncomp=ncomp)
      IF nuseind EQ 0 THEN message,'minx is less than minimum input x value'
      IF nuseind NE n_elements(xi) THEN BEGIN 
          npts = 100
          gauss[useind] = gauss[useind]/qgauss(gauss[useind], xi[useind], npts)
          gauss[comp] = 0.0
      ENDIF 
  ENDIF 

  return, gauss

end
