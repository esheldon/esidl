PRO binlumbynum_absrange, minAbsMag, maxAbsMag

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: binlumbynum_absrange, minAbsMag, maxAbsMag'
      return
  ENDIF 

  minAbsMag = [-22.0, -23.0, -24.0, -24.0, -25.0]
  maxAbsMag = [-12.0, -13.0, -14.0, -14.0, -14.0]

END 
