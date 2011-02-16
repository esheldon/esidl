PRO run_stripe_range_plots

    run_status = sdss_runstatus()
  w=where(run_status.stripe NE -1,nw)
  IF nw EQ 0 THEN message,'No good stripes'

  w2=rem_dup(run_status[w].stripe)
  w=w[w2]

  stripes = run_status[w].stripe
  nstripe=n_elements(stripes)
  FOR i=0L, nstripe-1 DO BEGIN
      stripe_range_plots, stripes[i], /dops
  ENDFOR 

END 
