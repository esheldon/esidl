FUNCTION wthetastruct, arrval

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: wtheta = wthetastruct(arrval)'
      return, -1
  ENDIF 

  ws = create_struct('nlenses', 0., $
                     'totpairs', 0., $
                     'binsize', 0., $
                     'rmin', 0., $
                     'rmax', 0., $
                     'h', 0., $
                     'meanr', arrval, $
                     'rmax_act', arrval, $
                     'npair', arrval, $
                     'tnpair', arrval,$
                     'area', arrval, $
                     'density', arrval)
  
  return, ws


END 
