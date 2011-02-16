PRO gauleg_idl, x1, x2, n, x, w

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: gauleg_idl, x1, x2, n, x, w'
      return
  ENDIF 

  x1 = double(x1)
  x2 = double(x2)
  n = long(n)
  x=dblarr(n)
  w=dblarr(n)

  EPS = double(3.e-11)

  m = (n+1)/2

  xm = (x1 + x2)/2d
  xl = (x2 - x1)/2d
  z1 = 0d

  FOR i=1, m DO BEGIN 
      z=cos( !dpi*(i-0.25)/(n+.5) )
      
      WHILE abs(z-z1) GT EPS DO BEGIN
          p1 = 1d
          p2 = 0d
          FOR j=1, n DO BEGIN
              p3 = p2
              p2 = p1
              p1 = ( (2d*j - 1d)*z*p2 - (j-1d)*p3)/j
          ENDFOR 
          pp = n*(z*p1 - p2)/(z^2 -1.)
          z1=z
          z=z1 - p1/pp
      ENDWHILE 
      x[i-1] = xm - xl*z
      x[n+1-i-1] = xm + xl*z
      w[i-1] = 2.*xl/( (1.-z^2)*pp*pp )
      w[n+1-i-1] = w[i-1]
  ENDFOR 

return
END 
      
