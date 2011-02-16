PRO check_admom, run, rerun, camcol, badfiles, front=front

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: check_admom, run, rerun, camcol, front=front'
      print
      print,' Default front is adatc'
      return
  ENDIF 


  fetch_dir, run, camcol, rerun, dir, corrdir=corrdir
  fetch_file_list,dir,files,fnums

  IF n_elements(front) EQ 0 THEN newfront = 'adatc' ELSE newfront=front

  rename_tsobj, files, corrdir, newfront, adatfiles, nchar

  fetch_file_list, corrdir, tfiles, front=newfront

  n_exp = n_elements(adatfiles)
  n_act = n_elements(tfiles)

  IF n_exp NE n_act THEN BEGIN
      print,'Number of files not right: ',n_exp,n_act
  ENDIF 

  i_act = 0
  delvarx, badfiles
  FOR i=0, n_exp-1 DO BEGIN

      IF adatfiles[i] NE adatfiles[i_act] THEN BEGIN
          print,'Bad: '+adatfiles[i]
          add_arrval, adatfiles[i], badfiles
      ENDIF ELSE BEGIN
          i_act = i_act+1
      ENDELSE 
  ENDFOR 

return 
END 
      



