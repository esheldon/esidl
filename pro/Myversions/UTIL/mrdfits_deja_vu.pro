FUNCTION mrdfits_deja_vu, file, hdr, hdr0=hdr0, silent=silent, deja_vu=deja_vu
  
  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: result=mrdfits_deja_vu(file, hdr, hdr0=hdr0, silent=silent, deja_vu=deja_vu'
      return,-1
  ENDIF 

  openr, lun, file, /get_lun, error=error
  IF error NE 0 THEN message,!ERR_STRING

  IF NOT keyword_set(deja_vu) THEN BEGIN 
      deja_vu=1
      tmp = mrdfits3(lun, 1, 0, hdr, silent=silent)
      free_lun, lun
      return,tmp
  ENDIF ELSE BEGIN 
      tmp = mrdfits3(lun, 1, 0, hdr, silent=silent, /deja_vu)
      free_lun, lun
      return,tmp
  ENDELSE 

END 
