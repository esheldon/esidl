FUNCTION denscont_probe, etan, ewidth,siginv,denscont

  ;; uses h=1.0, pc^2/Msolar
  model = 2.*denscont*siginv
  pofe = gaussprob(model, etan, ewidth)

  return, pofe

END 

PRO denscont_likelihood_mc, zlens, minzs, maxzs, npts, $
                            meanzs, sigzs, $
                            etan, etanerr, shapenoise, $
                            denscont, denscont_likelihood, $
                            silent=silent, $
                            prior_zs=prior_zs, prior_pofzs=prior_pofzs,$
                            justprior=justprior, meansigcritinv=meansigcritinv
                         
  ZZi = arrscl( dindgen(npts), minzs, maxzs )
  nrand = 400

  IF n_elements(prior_pofzs) NE 0 THEN BEGIN 

      IF keyword_set(justprior) THEN BEGIN 
          pofzsi = 1.
      ENDIF ELSE BEGIN 
          pofzsi = gaussprob(ZZi, meanzs, sigzs)
      ENDELSE 

      prior_pofzsi = interpol(prior_pofzs, prior_zs, ZZi)
      pofzsi = pofzsi*prior_pofzsi

      genrand, pofzsi, ZZi, nrand, randzs

  ENDIF ELSE BEGIN 
      randzs = randomu(seed, nrand, /normal)*sigzs + meanzs
  ENDELSE 

  siginv = sigmacritinv(zlens, randzs)

;  plothist, siginv/1.e-5,bin=0.1
;  key=get_kbrd(1)

  ewidth = sqrt( etanerr^2 + shapenoise^2 )

  ndenscont = n_elements(denscont)
  denscont_likelihood = dblarr(ndenscont)

  IF keyword_set(meansigcritinv) THEN BEGIN 
      mean_siginv = mean( siginv )
      denscont_likelihood[*] = denscont_probe(etan, ewidth, mean_siginv, denscont[*])
  ENDIF ELSE BEGIN 
      FOR i=0L, ndenscont-1 DO BEGIN 

          funcvals = denscont_probe(etan, ewidth, siginv, denscont[i])
          
          denscont_likelihood[i] = mean(funcvals)
          
      ENDFOR 
  ENDELSE 

END 
