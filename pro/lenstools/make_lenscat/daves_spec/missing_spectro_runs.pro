PRO missing_spectro_runs

  ;; print out which runs we are missing that have spectra.
  ;; based on the specto catalogs

  ;; note: when making these files we do matching by ra/dec
  ;; keeping those that did match to adatc files.  So this
  ;; may not be a perfect diagnostic.

  stripes = [9,10,11,12,30,31,32,33,34,35,36,37,42,43,76,82,86]

  nst = n_elements(stripes)

  w=where(!run_status.adatc_photo_v GE 5.3 AND $
          !run_status.adatc_photo_v LT 5.4)
  rmd = rem_dup(!run_status[w].run)
  rsruns = !run_status[w[rmd]].run

  FOR i=0L, nst-1 DO BEGIN 
      
      stripe = stripes[i]

      get_spectra_lcat,stripe,lcat,/spec

      runs = lcat[rem_dup(lcat.run)].run
      nr = n_elements(runs)
      match, runs, rsruns, m, mrs
      
      IF m[0] EQ -1 THEN BEGIN 
          nbad = n_elements(runs)
          nmiss = lonarr(nbad)
          FOR j=0L, nbad-1 DO BEGIN 
              wmiss = where(lcat.run EQ runs[j],tnmiss)
              nmiss[j] = tnmiss
          ENDFOR 
          print,'Stripe: '+stripe2string(stripe)
          print,'     Missing runs  #spectra'
          colprint,runs,nmiss
          print
      ENDIF ELSE IF n_elements(m) NE nr THEN BEGIN 
          remove, m, runs
          nbad = n_elements(runs)
          nmiss = lonarr(nbad)
          FOR j=0L, nbad-1 DO BEGIN 
              wmiss = where(lcat.run EQ runs[j],tnmiss)
              nmiss[j] = tnmiss
          ENDFOR 
          print,'Stripe: '+stripe2string(stripe)
          print,'     Missing runs  #spectra'
          colprint,runs,nmiss
          print
      ENDIF 

  ENDFOR 


END 
