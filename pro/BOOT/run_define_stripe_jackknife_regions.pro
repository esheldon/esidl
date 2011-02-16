PRO run_define_stripe_jackknife_regions, clusters=clusters

  IF keyword_set(clusters) THEN BEGIN 
      nPointsPerChunk = 100
      stripes = [ 9, 10, 11, 12, $
                 34, 35, 36, 37 ]
  ENDIF ELSE BEGIN 
      ;; not 16,26,42, 43, 44, 82 since don't pass current pixel masks
      nPointsPerChunk = 10
      stripes = [9,10,11,12,13,14,15,$
                 27,28,29,30,31,32,33,34,35,36,37,$
                 76,86]
  ENDELSE 


  nst = n_elements(stripes)

  FOR i=0L, nst-1 DO BEGIN
      stripe = stripes[i]
      define_stripe_jackknife_regions, stripe, nPointsPerChunk, $
        clusters=clusters
  ENDFOR 

END 
