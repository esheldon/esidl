PRO wthetarand_generate_etas, lambda, neta, angmax, eta, edgecut=edgecut

  COMMON wthetarand_block, seed, etarange, etarange_spectra, $
    mlambda, meta, FIRSTLAMBDA, LASTLAMBDA, Mlimit

  eta = double(randomu(seed, neta))
  maxiter = 100

  d2r = !dpi/180.0d0            ;Change degrees to radians.
  r2d = 180./!dpi               ;radians to degrees

  ;; etas generated from etarange not etarange_spectra
;  maxeta = interpol(etarange_spectra.maxeta, etarange_spectra.lambda, lambda)
;  mineta = interpol(etarange_spectra.mineta, etarange_spectra.lambda, lambda)
  maxeta = interpol(etarange.maxeta, etarange.lambda, lambda)
  mineta = interpol(etarange.mineta, etarange.lambda, lambda)


  totstr=ntostr(neta)
  FOR i=0L, neta-1 DO BEGIN 

      ;; IF i MOD 100L EQ 0 THEN print,ntostr(i+1)+'/'+totstr

      continue = 1
      teta = eta[i]
      iter = 0L
      WHILE continue DO BEGIN 
          teta = arrscl(teta, $
                        mineta[i], $
                        maxeta[i], $
                        arrmin=0., arrmax=1.)

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; make edgecut if requested
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          IF keyword_set(edgecut) THEN BEGIN 
              cosl = cos(lambda[i]*d2r)
              dmx = teta + angmax[i]/cosl
              dmn = teta - angmax[i]/cosl
              IF (dmx GT maxeta[i]) OR (dmn LT mineta[i]) THEN BEGIN 
                  teta = randomu(seed)
              ENDIF 
          ENDIF ELSE BEGIN 
              continue = 0
          ENDELSE   

;          gcirc, 0, teta*d2r, lambda[i]*d2r, $
;            meta, mlambda, dist
;          mdist = min(dist)*r2d
;          IF mdist GE angmax[i] THEN continue = 0 ELSE teta = randomu(seed)

          iter = iter+1
          IF iter GE maxiter THEN BEGIN
              teta = 1.e10
              continue = 0
          ENDIF 
      ENDWHILE 
      eta[i] = teta
  ENDFOR 

END 
