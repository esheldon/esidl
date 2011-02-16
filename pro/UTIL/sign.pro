FUNCTION sign, x

  nx = n_elements(x)

  IF nx EQ 1 THEN BEGIN 
      IF x[0] LT 0.0 THEN return, -1. ELSE return,1.
  ENDIF ELSE BEGIN 

      signs = replicate(1., nx)
      w=where(x LT 0.0,nw)
      IF nw NE 0 THEN signs[w] = -1.

      return,signs
  ENDELSE 

END 
