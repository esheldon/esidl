PRO zrandshear_random_points, stripe, nrand, wlens, angmax, $
                              llambda, mineta, maxeta, $
                              tsgals=tsgals, maskdir=maskdir, nomask=nomask,$
                              edgecut=edgecut

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: zrandshear_random_points, stripe, nrand, wlens, angmax, $'
      print,'         llambda, $'
      print,'         tsgals=tsgals, maskdir=maskdir, nomask=nomask'
      return
  ENDIF 

  COMMON rand_block, seed, $
    FIRSTLAMBDA, LASTLAMBDA, ETAMIN, ETAMAX, etarange, mask

  IF n_elements(etarange) EQ 0 THEN BEGIN 
      read_etarange, stripe, etarange, indir=etarangedir
  ENDIF 

  IF n_elements(mask) EQ 0 THEN BEGIN 
      read_stripe_mask, stripe, mask, tsgals=tsgals, indir=maskdir
  ENDIF 

  nstripe = n_elements(stripe)

  IF nstripe EQ 1 THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; these are special stripes. Spectroscopy is not a rectangle
      ;; in lambda, eta so to generate random points need etarange
      ;; and generate uniformly in lambda
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF stripe EQ 82 OR stripe EQ 10 OR stripe EQ 86 OR stripe EQ 76 THEN BEGIN 
          sinlambda = 0
      ENDIF ELSE BEGIN 
          sinlambda = 1
      ENDELSE 
  ENDIF ELSE BEGIN 
      message,'This program only works for n_elements(stripe) = 1'
  ENDELSE 

  sdssidl_setup,/silent
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; generate lambda's
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  ngood = 0L
  llambda = dblarr(nrand)
  
  bad = wlens
  nwlens = n_elements(wlens)
  nbad = nwlens
  print
  print,"Generating lambda's ",ntostr(nwlens)
  
  WHILE ngood LT nwlens DO BEGIN 
      
      IF ngood EQ 0 THEN print,'Generating ',ntostr(nbad) $
      ELSE print,'Re-Generating ',ntostr(nbad)
      
      ;; generate lambdas. 
      zrandshear_generate_lambdas, nbad, angmax[bad], tmplam, $
                                   sinlambda=sinlambda
      
      ;; Now we apply mask in rotated frame with new masks
      
      ;; now apply mask
      IF NOT keyword_set(nomask) THEN BEGIN 
          apply_mask, mask, tmplam, tbad, tgood, edgecut=angmax[bad]
      
          IF tgood[0] NE -1 THEN BEGIN 
              ;; we found some good ones
              good = bad[tgood]
              IF tbad[0] NE -1 THEN BEGIN 
                  bad = bad[tbad]
                  nbad = n_elements(bad)
              ENDIF ELSE nbad=0L
              llambda[good] = tmplam[tgood]
              
              ntgood = n_elements(tgood)
              ngood = ngood + ntgood
          ENDIF 
      ENDIF ELSE BEGIN
          llambda = tmplam
          ngood = nwlens
      ENDELSE 
  ENDWHILE 
  ;; lambda's emerge unsorted and rotated
  
  ;; sort the subscript array
  s = sort(llambda[wlens])
  wlens = wlens[s]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; these are special stripes. Spectroscopy is not a rectangle
  ;; in lambda, eta so to generate random points need etarange
  ;; and generate uniformly in lambda
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT sinlambda THEN BEGIN 
      
      ;; generate range for eta
      maxeta = interpol(etarange.maxeta, etarange.clambda, llambda)
      mineta = interpol(etarange.mineta, etarange.clambda, llambda)

      ;; alter maxeta,mineta for edge cuts
      IF keyword_set(edgecut) THEN BEGIN 

          maxeta = maxeta - angmax
          mineta = mineta + angmax

          ;; Did user input redshifts that are too small?
          wbad = where(mineta GE maxeta, nbad)
          IF nbad NE 0 THEN message,'Problem with mineta, maxeta'

      ENDIF 

  ENDIF ELSE BEGIN 

      primary_bound, stripe, bound

      maxeta = replicate(bound.etamax, nrand)
      mineta = replicate(bound.etamin, nrand)

      IF keyword_set(edgecut) THEN BEGIN 

          maxeta2 = interpol(etarange.maxeta, etarange.clambda, llambda) - angmax
          mineta2 = interpol(etarange.mineta, etarange.clambda, llambda) + angmax

          maxeta = maxeta2 < maxeta
          mineta = mineta2 > mineta
          
          ;; Did user input redshifts that are too small?
          wbad = where(mineta GE maxeta, nbad)
          IF nbad NE 0 THEN message,'Problem with mineta, maxeta'

      ENDIF 

  ENDELSE 




END 
