PRO rotate_ra, ra

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: rotate_ra, ra'
      return
  ENDIF 

  ra = ra + 180d
  w=where(ra GT 360d, nw)
  IF nw NE 0 THEN ra[w] = ra[w] - 360d

  return
END
