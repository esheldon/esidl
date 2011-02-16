;; This is ours, a simplified unique identifier

FUNCTION sdss_spectroid, plate, mjd, fiberid

  ten=ulong64(10)

  p1 = 0L
  p2 = 4L
  p3 = 10L

  np = n_params()
  IF np EQ 1 THEN BEGIN
      super = $
        ulong64(plate.fiberid > 0)*ten^p1 + $
        ulong64(plate.mjd > 0)*ten^p2 + $
        ulong64(plate.plate > 0)*ten^p3
  ENDIF ELSE IF np EQ 3 THEN BEGIN 
      
      np = n_elements(plate)
      nj = n_elements(mjd)
      nf = n_elements(fiberid)

      IF total_int([np, nj, nf]) NE 3*np THEN BEGIN 
          print,'All arrays must be same size'
          return,-1
      ENDIF 

      super = $
        ulong64(fiberid > 0)*ten^p1 + $
        ulong64(mjd > 0)*ten^p2 + $
        ulong64(plate > 0)*ten^p3
  ENDIF ELSE BEGIN 
      print,'-Syntax: superid=sdss_specid(plate, mjd, fiberid -OR- struct)'
      return,-1
  ENDELSE 
  return,super

END 
