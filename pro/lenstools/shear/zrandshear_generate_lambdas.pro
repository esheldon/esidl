PRO zrandshear_generate_lambdas, nlam, angmax, lambda, sinlambda=sinlambda

  COMMON rand_block, seed, $
    FIRSTLAMBDA, LASTLAMBDA, ETAMIN, ETAMAX, etarange, mask, use_etarange

  IF keyword_set(sinlambda) THEN BEGIN 

      ;; generate in sin(lam): for when using primary bounds

      sinflam = sin( (FIRSTLAMBDA + angmax)*!d2r )
      sinllam = sin( (LASTLAMBDA  - angmax)*!d2r )

      slambda = randomu(seed, nlam, /double)

      FOR i=0L, nlam-1 DO BEGIN 
          
          slambda[i] = arrscl(slambda[i], $
                              sinflam[i], $
                              sinllam[i], $
                              arrmin=0., arrmax=1.)
      ENDFOR 

      lambda = asin(slambda)*!r2d

  ENDIF ELSE BEGIN 
      lambda = randomu(seed,nlam, /double)

      FOR i=0L, nlam-1 DO BEGIN 
          
          lambda[i] = arrscl(lambda[i], $
                             FIRSTLAMBDA + angmax[i], $
                             LASTLAMBDA  - angmax[i], $
                             arrmin=0., arrmax=1.)
      ENDFOR 
  ENDELSE 

END 
