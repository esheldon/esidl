PRO denscont_likelihood, zlens, minzs, maxzs, npts, $
                         meanzs, sigzs, $
                         etan, etanerr, shapenoise, $
                         denscont, denscont_likelihood, $
                         mdenscont_likelihood, $
                         mean_denscont, err_denscont, $
                         silent=silent, $
                         prior_zs=prior_zs, prior_pofzs=prior_pofzs,$
                         justprior=justprior, meansigcritinv=meansigcritinv
                         
  gauleg, minzs, maxzs, npts, ZZi, WWi

  IF n_elements(prior_pofzs) NE 0 THEN BEGIN 

      IF keyword_set(justprior) THEN BEGIN 
          pofzsi = 1.
      ENDIF ELSE BEGIN 
          pofzsi = gaussprob(ZZi, meanzs, sigzs)
      ENDELSE 

      prior_pofzsi = interpol(prior_pofzs, prior_zs, ZZi)
      pofzsi = pofzsi*prior_pofzsi

  ENDIF ELSE BEGIN 
      pofzsi = gaussprob(ZZi, meanzs, sigzs)
  ENDELSE 

  ;; renormalize
  norm = total(pofzsi*WWi)
  pofzsi = pofzsi/norm

  siginv = sigmacritinv(zlens, ZZi)

  ewidth = sqrt( etanerr^2 + shapenoise^2 )

  ndenscont = n_elements(denscont)
  denscont_likelihood = dblarr(ndenscont)
  mdenscont_likelihood = denscont_likelihood

  ;; mean
  func = pofzsi*siginv
  mean_siginv = total(WWi*func)

  mean_denscont = etan/2./mean_siginv
  err_denscont = ewidth/2./mean_siginv

  mdenscont_likelihood[*] = denscont_probe(etan, ewidth, mean_siginv, denscont[*])

;  varfunc = pofzsi*(siginv-mean_siginv)^2
;  siginv_var = total( WWi*varfunc )
  
;  ewidth_send = sqrt(ewidth^2 + denscont[*]^2*siginv_var)
;  mdenscont_likelihood[*] = denscont_probe(etan, ewidth_send, mean_siginv, denscont[*])

  IF keyword_set(meansigcritinv) THEN BEGIN 
      ;; just do the mean
      denscont_likelihood[*] = mdenscont_likelihood[*]
      return
  ENDIF ELSE BEGIN 
      ;; do full calculation
      FOR i=0L, ndenscont-1 DO BEGIN 
          
          funcvals = pofzsi*denscont_probe(etan, ewidth, siginv, denscont[i])
          
          denscont_likelihood[i] = total(funcvals*WWi)
      
      ENDFOR 
  ENDELSE 

END 
