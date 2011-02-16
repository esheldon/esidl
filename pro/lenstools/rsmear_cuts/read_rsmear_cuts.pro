PRO read_rsmear_cuts, run, rerun, purity, clr, rcutstr, status=status, camcol=camcol, hirata=hirata

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: read_rsmear_cuts, run, rerun, purity, clr, rcutstr'
      return
  ENDIF 

  rsmear_cuts_files, run, rerun, purity, clr, fitfile, $
                     camcol=camcol, hirata=hirata
  
  print
  print,'Reading rsmear_cuts file: ',fitfile
  rcutstr = mrdfits(fitfile, 1, /silent, status=status)

  IF status NE 0 THEN print,'No such rsmear cuts file: ',fitfile

END 
