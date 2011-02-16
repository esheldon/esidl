PRO vagc_run_setuplens, random=random

;  rmin = 20.0
;  rmax = 10000.0
;  nbin_OR_binsize = 18
;  logbin = 1

  rmin = 20.0
  rmax = 1020.0
  nbin_OR_binsize = 100.
  logbin = 0

  ;; No point in 82..
  stripes = [9,10,11,12,13,14,15, $
             27,28,29,30,31,32,33,34,35,36,37, $
             76, 86]

  compcut = 0.9

  IF NOT keyword_set(random) THEN BEGIN 

      vagc_setuplens, rmin, rmax, nbin_OR_binsize, stripes, $
        compcut=compcut, logbin=logbin

  ENDIF ELSE BEGIN 
      
      randnum = lindgen(10)
      nrand = n_elements(randnum)

      vagc_setuprand, randnum, rmin, rmax, nbin_OR_binsize, stripes, $
        compcut=compcut, logbin=logbin

  ENDELSE 

END 
