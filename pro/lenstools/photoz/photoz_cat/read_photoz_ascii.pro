;; For use with Budavari photozs

PRO read_photoz_ascii, file, struct, status=status

  status=1
  IF n_params() LT 1 THEN BEGIN  
      print,'-Syntax: read_photoz_ascii, file, struct'
      return
  ENDIF 

  IF NOT fexist(file) THEN BEGIN
      print,'No such file: '+file
      return
  ENDIF 

  print,'Reading photoz file: '+file

  nlines = numlines(file)

  make_photoz_struct, tstr

  struct = replicate(tstr, nlines)

  openr, lun, file, /get_lun
  readf, lun, struct
  free_lun, lun

  status=0

END 
