
FUNCTION mean_sigmacritinv, zlens, meanzs, sigzs, npts, $
                            prior_zs=prior_zs,  prior_pofzs=prior_pofzs

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: result = mean_sigmacritinv(zlens, meanzs, sigzs, npts, $'
      print,'                           prior_zs=prior_zs,  prior_pofzs=prior_pofzs)'
      print,'/use_lambda is default'
      return,-1.
  ENDIF 

  nsig = 4.0
  minzs = meanzs - nsig*sigzs
  maxzs = meanzs + nsig*sigzs

  gauleg, minzs, maxzs, npts, ZZi, WWi

  IF n_elements(prior_pofzs) NE 0 THEN BEGIN 

      pofzsi = gaussprob(ZZi, meanzs, sigzs)

      prior_pofzsi = interpol(prior_pofzs, prior_zs, ZZi)
      pofzsi = pofzsi*prior_pofzsi

      ;; renormalize
      norm = total(pofzsi*WWi)
      pofzsi = pofzsi/norm

  ENDIF ELSE BEGIN 
      ;; already normalized
      pofzsi = gaussprob(ZZi, meanzs, sigzs)
  ENDELSE 

  siginv = sigmacritinv(zlens, ZZi)
                   

  func = pofzsi*siginv
  mean_siginv = total(WWi*func)

  return, mean_siginv


END 
