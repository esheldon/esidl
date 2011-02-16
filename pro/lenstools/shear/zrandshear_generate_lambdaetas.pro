PRO zrandshear_generate_lambdaetas, num, lambda, eta

  COMMON mask_block, seed, $
    FIRSTLAMBDA, LASTLAMBDA, ETAMIN, ETAMAX, $
    etarange, mask, use_etarange_edgecut, special

  IF special THEN BEGIN 

      ;; generate uniformly in lambda, eta since using whole width
      ;; of scan
  
      IF float(!version.release) GE 5.5 THEN BEGIN 
          lambda = arrscl( randomu(seed, num, /double), $
                           FIRSTLAMBDA, LASTLAMBDA, $
                           arrmin=0., arrmax=1.)
          
          eta = randomu(seed, num, /double)
      ENDIF ELSE BEGIN 
          lambda = dblarr(num) & eta = lambda
          
          lambda[*] = arrscl( randomu(seed, num), $
                           FIRSTLAMBDA, LASTLAMBDA, $
                           arrmin=0., arrmax=1.)
          
          eta[*] = randomu(seed, num)
      ENDELSE 
      maxeta = interpol(etarange.maxeta, etarange.clambda, lambda)
      mineta = interpol(etarange.mineta, etarange.clambda, lambda)

      FOR i=0L, num-1 DO BEGIN 
          eta[i] = arrscl(eta[i], $
                          maxeta[i], mineta[i], $
                          arrmin=0., arrmax=1.)
      ENDFOR 

  ENDIF ELSE BEGIN 
      
      ;; generate uniformly in sin(lambda) and eta
      sinflam = sin( FIRSTLAMBDA*!d2r )
      sinllam = sin( LASTLAMBDA*!d2r )

      IF float(!version.release) GE 5.5 THEN BEGIN 
          slambda = arrscl( randomu(seed, num, /double), sinflam, sinllam, $
                            arrmin=0., arrmax=1.)
      
          lambda = asin(temporary(slambda))*!r2d

          eta = arrscl(randomu(seed, num, /double), $
                       ETAMIN, ETAMAX, $
                       arrmin=0., arrmax=1.)
      ENDIF ELSE BEGIN 
          lambda = dblarr(num) & eta = lambda


          slambda = double(arrscl( randomu(seed, num), sinflam, sinllam, $
                                   arrmin=0., arrmax=1.) )
      
          lambda[*] = asin(temporary(slambda))*!r2d

          eta[*] = arrscl(randomu(seed, num), $
                          ETAMIN, ETAMAX, $
                          arrmin=0., arrmax=1.)
      ENDELSE 
  ENDELSE 

END 
