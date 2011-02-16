PRO admomascii2fits, infile, outfile, status=status

  ;; status is 1 unless we reach end
  status=1

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: admomascii2fits, infile [, outfile, status=status]'
      return
  ENDIF 

  outfile = repstr(infile,'.dat','.fit')

  print
  print,'Infile = ',infile
  print,'Outfile = ',outfile
  print

  read_admom, infile, struct

  mwrfits, struct, outfile, /create

  status=0

END 
