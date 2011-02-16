PRO run_make_scat53, runs, reruns, clr, hirata=hirata

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: run_make_scat53, runs, reruns, clr, hirata=hirata'
      return
  ENDIF 

  nrun = n_elements(runs)

  FOR i=0L, nrun-1 DO BEGIN 
      
      make_scat53_meane, runs[i], reruns[i], clr, hirata=hirata
          
  ENDFOR 

END 
