PRO testsurvey2eq, scat
 
  t=systime(1)
  survey2eq_old, scat.lambda, scat.eta, ra1, dec1
  ptime,systime(1)-t

  t=systime(1)
  survey2eq, scat.lambda, scat.eta, ra2, dec2
  ptime,systime(1)-t

  print,max( abs(ra2-ra1) ),max( abs(ra2-ra1) )

  return
END 
