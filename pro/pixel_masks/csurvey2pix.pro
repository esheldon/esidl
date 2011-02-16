FUNCTION csurvey2pix, clambda, ceta, resolution=resolution

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: pixnum = csurvey2pix(clambda, ceta, resolution=)'
      print
      print,"just a wrapper for Ryan's ang2pix"
      return,-1
  ENDIF 
  ang2pix, clambda, ceta, index_x, index_y, pixnum, /survey, $
    resolution=resolution

  return,pixnum

END 
