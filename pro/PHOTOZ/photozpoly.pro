FUNCTION photozpoly, u, g, r, i

  IF n_params() LT 4 THEN BEGIN
      print,'-Syntax: result=photozpoly(u,g,r,i)'
      print,' error ~ 0.035 to z=.25'
      return,0.
  ENDIF 

  rr = r^2
  gg = g^2
  uu = u^2
  ii = i^2
  rg = r*g
  ru = r*u
  gu = g*u
  gi = g*i
  ui = u*i

  return, -2.465 - 5.603*r - 1.095*g + 0.9343*u + 6.049*i $
                 + 0.5981*rr + 0.3944*rg - 0.4927*gg - 0.8412*ru $
                 + 0.5700*gu - 0.06564*uu - 0.3275*r*i $
                 + 0.07044*gi + 0.3421*ui - 0.2563*ii

END 
