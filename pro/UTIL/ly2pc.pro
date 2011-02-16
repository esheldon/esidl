PRO ly2pc, ly, pc

  IF n_params() LT 1 THEN BEGIN
      print,'Syntax: ly2pc, lydist, pcdist
      print,' Give dist in light years, return dist in parcecs'
      return
  ENDIF 

  pc = ly/3.26

  return
END 
