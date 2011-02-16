FUNCTION wthetalumwsumstruct, arrval,beta

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: wtheta = wthetasumstruct(arrval,betaarray)'
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
                     'beta',beta, $
                     'rsum', arrval, $
                     'lsum',arrval,$
                     'lerrsum1',arrval,$
                     'lerrsum2',arrval,$
                     'lerrsum3',arrval,$
                     'lwsum',arrval,$
                     'rmax_act', arrval, $
                     'rmin_act', arrval, $
                     'npsum', arrval2, $
                     'wsum', arrval)
  
  return, ws


END 
