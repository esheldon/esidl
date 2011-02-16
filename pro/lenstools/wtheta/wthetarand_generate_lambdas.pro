PRO wthetarand_generate_lambdas, nlam, angmax, lambda

  COMMON wthetarand_block, seed, etarange, etarange_spectra, $
    mlambda, meta, FIRSTLAMBDA, LASTLAMBDA, Mlimit

  lambda = double(randomu(seed,nlam))

  FOR i=0L, nlam-1 DO BEGIN 

      lambda[i] = arrscl(lambda[i], $
                         FIRSTLAMBDA + angmax[i], $
                         LASTLAMBDA  - angmax[i], $
                         arrmin=0., arrmax=1.)
  ENDFOR 

END 
