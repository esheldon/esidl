FUNCTION areapoly, r, p

  f = r
  f[*] = 1.0

  np = n_elements(p)
  FOR i=0,np-1 DO f[*] = f[*] + p[i]*r[*]^(i+1.0)

  return, f

END 
