FUNCTION ufunc, Rmax, l=l, cen=cen
  
  IF n_params() EQ 0 THEN BEGIN
      print,'-Syntax: result = ufunc(Rmax, l=l, cen=cen)'
      return,-1
  ENDIF 

  Rmax2 = long(Rmax)^2
  size = long(2.*Rmax)+1
  u = fltarr(size, size)
  
  IF n_elements(l) EQ 0 THEN l=1.
  IF n_elements(cen) EQ 0 THEN cen = [(size-1.)/2.,(size-1.)/2.]

  cx = cen(0)
  cy = cen(1)

  index = lindgen(size^2)

  x = index MOD size
  y = index/size
  
  r2 =  (x-cx)^2 + (y-cy)^2

  w=where( r2 LE Rmax2, nw)
  IF nw EQ 0 THEN BEGIN
      print,'No good'
      return, -1
  ENDIF 
  fac = (1./!pi/Rmax2)*(l+2.)^2
  u[w] = fac*(1. - r2[w]/Rmax2)^l*(1./(l+2.) - r2[w]/Rmax2)
  return, u

END 
