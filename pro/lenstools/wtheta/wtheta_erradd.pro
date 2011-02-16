PRO wtheta_erradd, file, lensumfile

  print
  print,'Input file: ',file
  print,'Input lensum file: ',lensumfile
  print

  l=mrdfits(file,1)
  lensum=mrdfits(lensumfile,1)

  wtheta_calcmeanerr, lensum, meanlumerr

  add_tag, l, 'meanlumerr', meanlumerr, newl

  ;; is it a random file or regular file?

  tmp = str_sep(file, 'wthetalumw')
  IF n_elements(tmp) EQ 1 THEN BEGIN 
      ;; was random file
      newfile=repstr(file, 'wthetarandlumw', 'wthetarandlumwerradd')
  ENDIF ELSE BEGIN 
      ;; was regular file
      newfile = tmp[0]+'wthetalumwerradd'+tmp[1]
  ENDELSE 

  print
  print,'Input file: ',file
  print,'Input lensum file: ',lensumfile
  print,'New output file: ',newfile
  print

  mwrfits, newl, newfile, /create

END 
