FUNCTION run2stripenum, run

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: stripe = run2stripenum(run)'
      print,' run can be an array'
      return,-1
  ENDIF 

  COMMON runlist_block, runstruct

  IF n_elements(runstruct) EQ 0 THEN read_runlist, runstruct

  nrun = n_elements(run)
  IF nrun EQ 1 THEN BEGIN 
      wrun = where(runstruct.run EQ run, nwrun)
      IF nwrun EQ 0 THEN return,-1 ELSE return,runstruct[wrun].stripe
  ENDIF ELSE BEGIN 
      stripes = replicate(-1, nrun)
      match, run, runstruct.run, mrun, mstruct
      IF mstruct[0] NE -1 THEN BEGIN 
          stripes[mrun] = runstruct[mstruct].stripe
      ENDIF 
      return,stripes
  ENDELSE 
  

END 
