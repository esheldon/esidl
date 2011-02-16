PRO testeq2survey, scat

  t=systime(1)
  eq2survey_old, scat.ra, scat.dec, lam1, eta1
  ptime,systime(1)-t

  t=systime(1)
  eq2survey, scat.ra, scat.dec, lam2, eta2
  ptime,systime(1)-t

  print,max( abs(lam2-lam1) ),max( abs(eta2-eta1) )

  return
END 
