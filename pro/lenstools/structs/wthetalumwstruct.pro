FUNCTION wthetalumwstruct, arrval, beta

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: wtheta = wthetastruct(arrval,betaarray)'
      return, -1
  ENDIF 

  nbeta=n_elements(beta)
  arrval2 = fltarr(n_elements(arrval), nbeta)

  ws = create_struct('nlenses', 0., $
                     'totpairs', 0., $
                     'binsize', 0., $
                     'rmin', 0., $
                     'rmax', 0., $
                     'h', 0., $
                     'beta', beta, $
                     'meanr', arrval, $
                     'meanlum',arrval,$
                     'meanlumerr',arrval,$
                     'tmeanlum', arrval,$
                     'tmeanlumerr', arrval,$
                     'rmax_act', arrval, $
                     'rmin_act', arrval, $
                     'npair', arrval2, $
                     'tnpair', arrval2, $
                     'area', arrval, $
                     'density', arrval2)
  
  return, ws


END 
