PRO randomize_lameta,t

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: randomize_lameta,t'
      return
  ENDIF 

  print,'Randomizing lambdas'

  COMMON wthetarand_block, seed, etarange, etarange_spectra, $
    mlambda, meta, FIRSTLAMBDA, LASTLAMBDA, Mlimit

  seed=long(systime(1))
  read_etarange, 10, etarange

  nn=n_elements(t)

  IF NOT tag_exist(t,'lambda') THEN BEGIN 
      print,'adding lambda/eta tags'
      add_tags, t, ['lambda','eta'], ['0d','0d'], newt
      delvarx, t
      t=temporary(newt)

      eq2survey, t.ra, t.dec, llambda, leta
      t.lambda = llambda
      t.eta = leta
  ENDIF 

  maxlam=max(t.lambda, min=minlam)

  lambda = double( arrscl(randomu(seed, nn), minlam, maxlam, $
                          arrmin=0., arrmax=1.) )

  wthetarand_generate_etas, lambda, nn, fltarr(nn), eta

  print,'Sorting lambdas'
  s=sort(lambda)

  t.lambda = lambda
  t.eta = eta

  t=temporary(t[s])

;  plot,t.lambda,t.eta,psym=3

END 
