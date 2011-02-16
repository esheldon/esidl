FUNCTION arange, array, $
                 subscript_min=subscript_min, $
                 subscript_max=subscript_max, $
                 NaN=NaN, $
                 dimension=dimension

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: rng = arange(array, $'
      print,'                      subscript_min=subscript_min, $'
      print,'                      subscript_max=subscript_max, $'
      print,'                      NaN=NaN, $'
      print,'                      dimension=dimension)'
      return,-1
  ENDIF 

  amin=min(array, subscript_min, $
           max=amax, $
           subscript_max=subscript_max, NaN=NaN, $
           dimension=dimension)
  return, [amin, amax]

END 
