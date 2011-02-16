FUNCTION qromb_pzs, zs

  COMMON dqromb_blk, zL, meanzsi, sigzsi, priorzsi, priorpofzsi, etani, ewidthi, tdenscont
  return, interpol(priorpofzsi, priorzsi, zs)*gaussprob(zs, meanzsi, sigzsi)

END 

FUNCTION qromb_scinv, zs

  COMMON dqromb_blk, zL, meanzsi, sigzsi, priorzsi, priorpofzsi, etani, ewidthi, tdenscont

  pzs = qromb_pzs(zs)

  return, sigmacritinv(zl, zs)*pzs

END 

FUNCTION qromb_like, zs

  COMMON dqromb_blk, zL, meanzsi, sigzsi, priorzsi, priorpofzsi, etani, ewidthi, tdenscont

  siginv = sigmacritinv(zl, zs)
  return, qromb_pzs(zs)*denscont_probe(etani, ewidthi, siginv, tdenscont)

END 

PRO denscont_likelihood_qromb, zlens, minzs, maxzs, npts, $
                               meanzs, sigzs, $
                               etan, etanerr, shapenoise, $
                               denscont, denscont_likelihood, $
                               silent=silent, $
                               prior_zs=prior_zs, prior_pofzs=prior_pofzs,$
                               justprior=justprior, meansigcritinv=meansigcritinv
                         

  COMMON dqromb_blk, zL, meanzsi, sigzsi, priorzsi, priorpofzsi, etani, ewidthi, tdenscont

  zL = zlens
  meanzsi = meanzs
  sigzsi = sigzs
  priorzsi = prior_zs
  priorpofzsi = prior_pofzs

  eps = 1.e-5

  wset,3
  print,meanzs
  tzsi = arrscl( findgen(npts), minzs, maxzs )
  tscinv = qromb_scinv(tzsi)
  tpofzsi = qromb_pzs(tzsi)
  plot,tzsi,tpofzsi
  plot,tzsi,tscinv
  wset,0

  ;; Get norm of pofzsi tot

  norm = qromb('qromb_pzs', minzs, maxzs, eps=eps)

  ewidth = sqrt( etanerr^2 + shapenoise^2 )
  ewidthi = ewidth
  etani = etan

  ndenscont = n_elements(denscont)
  denscont_likelihood = dblarr(ndenscont)

  IF keyword_set(meansigcritinv) THEN BEGIN 
      mean_siginv = qromb('qromb_scinv', minzs, maxzs, eps=eps)/norm
      denscont_likelihood[*] = denscont_probe(etan, ewidth, mean_siginv, denscont[*])
  ENDIF ELSE BEGIN 

      FOR i=0L, ndenscont-1 DO BEGIN 
          
          tdenscont = denscont[i]
          print,'i=',i
          denscont_likelihood[i] = qromb('qromb_like', minzs, maxzs, eps=eps)
      
      ENDFOR 
      denscont_likelihood = denscont_likelihood/norm
  ENDELSE 

END 
