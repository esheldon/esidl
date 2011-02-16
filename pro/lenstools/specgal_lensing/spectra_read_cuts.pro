PRO read_cuts, range_struct

  IF n_params() LT 1 THEN BEGIN 
      print,'read_cuts, range_struct'
  ENDIF 

  file = '~/lensout/cuts.fit'
  print,'Reading file: ',file

  range_struct = mrdfits(file, 1)

END 
