PRO read_selection_func, struct

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: read_selection_func, struct'
      return
  ENDIF 

  fitfile = '~/lensout/selection_function.fit'
  
  IF NOT fexist(fitfile) THEN BEGIN 
      file = '~/lensout/selection_function.dat'
      
      readcol, file, cz, phi
      
      struct = create_struct('cz', cz, $
                             'z', cz/3.e5, $
                             'phi', phi)

      mwrfits, struct, fitfile, /create
  ENDIF ELSE struct = mrdfits(fitfile, 1)


END 
