PRO zrandshear_random_points_sphpix_generate, num, bound_array, lambda, eta

  ;; generate uniformly in sin(lambda) and eta
  sinflam = sin( bound_array.lammin*!d2r )
  sinllam = sin( bound_array.lammax*!d2r )
  
  IF float(!version.release) GE 5.5 THEN BEGIN 
      slambda = arrscl( randomu(seed, num, /double), $
                        sinflam, sinllam, $
                        arrmin=0., arrmax=1.)
      
      lambda = asin(temporary(slambda))*!r2d
      
      eta = arrscl(randomu(seed, num, /double), $
                   bound_array.etamin, bound_array.etamax, $
                   arrmin=0., arrmax=1.)

  ENDIF ELSE BEGIN 
      lambda = dblarr(num) & eta = lambda
      
      
      slambda = double(arrscl( randomu(seed, num), sinflam, sinllam, $
                               arrmin=0., arrmax=1.) )
      
      lambda[*] = asin(temporary(slambda))*!r2d
      
      eta[*] = arrscl(randomu(seed, num), $
                      bound_array.etamin, bound_array.etamax, $
                      arrmin=0., arrmax=1.)
  ENDELSE 

END 

PRO zrandshear_random_points_sphpix, stripes, nrand, angmax, llambda, leta, $
                                      edgeflag, silent=silent, $
                                      noetacut=noetacut, etaor=etaor

  ;; angmax,edgeflag unused until I implement the edge checking

  primary_bound_multi, stripes, bound_array, overall_bound

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
                      
      IF nbad EQ 0 THEN stop
      zrandshear_random_points_sphpix_generate, nbad, overall_bound, $
                                                tmplam,tmpeta

      ;; TEMPORARY EDGEFLAG STUFF
      tedgeflag = tmplam

      ;; First apply spectroscopic mask
      tbadlist = bytarr(nbad)
      apply_sphpoly_mask, tmplam, tmpeta, tbad, tgood1, /lameta, compcut=0
      IF tbad[0] NE -1 THEN tbadlist[tbad] = 1b

      IF tgood1[0] NE -1 THEN BEGIN 

          ;; now apply the pixel mask (just boundary for now)
          apply_pixel_mask, tmplam[tgood1], tmpeta[tgood1], tbad2, tgood2

          ;; TO DO: write edge-checking routine in C
          ;; for the lensing stuff (not needed for Amy's stuff)

          IF tgood2[0] NE -1 THEN BEGIN 
          
              tgood = tgood1[tgood2]
              IF tbad2[0] NE -1 THEN tbadlist[tgood1[tbad2]] = 1b
              tbad = where(tbadlist EQ 1b)

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
      ENDIF 
  ENDWHILE 

  ;; lambda's emerge unsorted



END 
