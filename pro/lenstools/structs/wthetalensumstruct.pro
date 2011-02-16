FUNCTION wthetalensumstruct, arrval

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: wtheta = wthetasumstruct(arrval)'
      return, -1
  ENDIF

  ws = create_struct('zindex', 0L, $
                     'totpairs', 0., $
                     'binsize', 0., $
                     'rmin', 0., $
                     'rmax', 0., $
                     'h', 0., $
                     'rsum', arrval, $
                     'rmax_act', arrval, $
                     'npsum', arrval,$
                     'wsum',arrval)
  
  return, ws


END 
