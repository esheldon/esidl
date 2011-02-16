PRO make_process_reruns, runs, reruns, stripes, bad

  ;; spit out runs that should be processed, but are not yet

  IF n_params() EQ 0 THEN BEGIN 
      print,'-Syntax: make_process_reruns, runs, reruns, stripes, bad'
      return
  ENDIF 

  delvarx,runs,reruns,stripes,bad

  ;; taken from /sdss4/data1/esheldon/GAL_GAL/spectra/fixcollated_12_15_00.fit
  cruns=[ 94, $
         125, $
         259, $
         273, $
         745, $
         752, $
         756, $
        1033, $
        1035, $
        1043, $
        1056, $
        1140, $
        1231, $
        1331, $
        1336, $
        1339, $
        1345, $
        1350, $
        1356, $
        1359, $
        1402, $
        1450]

  doneruns = [1140, $
              1231, $
              1241, $
              1331, $
              1345, $
              752, $
              756, $
              745, $
              1350, $
              1336, $
              1339, $
              1356, $
              1359, $
              1402, $
              1450, $
              94, $
              125]


  sdssidl_setup,/silent

  w=where(!run_status.bad EQ 0)
  trunstat = !run_status[w]

  nrun = n_elements(cruns)
  
  FOR i=0L, nrun-1 DO BEGIN 

      w=where(trunstat.run EQ cruns[i], nw)
      w2 = where(doneruns EQ cruns[i], nw2)

      IF (nw GT 0) AND (nw2 EQ 0) THEN BEGIN 
          add_arrval, cruns[i], runs
          add_arrval, max( trunstat[w].rerun ), reruns
          add_arrval, trunstat[w[0]].stripe, stripes
          add_arrval, trunstat[w[0]].bad, bad
      ENDIF 
  ENDFOR 
              
  colprint,runs,reruns,stripes,bad


END 
