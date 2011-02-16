PRO zrandshear_random_points, nrand, angmax, llambda, leta, edgeflag, silent=silent, etaor=etaor, noetacut=noetacut

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: zrandshear_random_points, nrand, angmax, llambda, leta, edgeflag, silent=silent, etaor=etaor, noetacut=noetacut'
      return
  ENDIF 

  COMMON mask_block, seed, $
    FIRSTLAMBDA, LASTLAMBDA, ETAMIN, ETAMAX, etarange, mask, use_etarange

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; generate clambda's and ceta's
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(etaor) THEN print,'ORing the eta checks'
  IF keyword_set(noetacut) THEN print,'No eta checks at all'

  ngood = 0L

  llambda = dblarr(nrand)
  leta    = dblarr(nrand)
  edgeflag = bytarr(nrand)

  bad = lindgen(nrand)

  nbad = nrand
  IF NOT keyword_set(silent) THEN BEGIN 
      print
      print,"Generating lambda's and eta's: ",ntostr(nrand)
  ENDIF 
  WHILE ngood LT nrand DO BEGIN 
      
      IF NOT keyword_set(silent) THEN BEGIN 
          IF ngood EQ 0 THEN print,'Generating ',ntostr(nbad) $
          ELSE print,'Re-Generating ',ntostr(nbad)
      ENDIF 
      ;; generate lambdas. 
      zrandshear_generate_lambdaetas, nbad, tmplam, tmpeta
      
      ;; now apply mask
      zobjshear_apply_mask, tmplam, tmpeta, angmax[bad], $
                            tbad, tgood, tedgeflag, $
                            etaor=etaor, noetacut=noetacut

      IF tgood[0] NE -1 THEN BEGIN 
          
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; we found some good ones: this will subscript the
          ;; new good ones in llambda, leta
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          good_ids = bad[tgood]

          ;; copy in the good ones
          llambda[good_ids] = tmplam[tgood]
          leta[good_ids]    = tmpeta[tgood]
          edgeflag[good_ids] = tedgeflag[tgood]

          ;; how many good ones now?
          ngood = ngood + n_elements(tgood)

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; now set bad ones to the ones that failed this
          ;; time around, or to zero if we are finished
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF tbad[0] NE -1 THEN BEGIN 
              bad = bad[tbad]
              nbad = n_elements(bad)
          ENDIF ELSE nbad=0L
          
      ENDIF 

  ENDWHILE 

  ;; lambda's emerge unsorted
  
END 
