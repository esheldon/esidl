PRO maxsis_s2n, zlens, zsource, SN0=SN0, nlenses=nlenses, oplot=oplot, $
                sigrange=sigrange, yrange=yrange

  IF n_elements(zlens) EQ 0 THEN zlens=.15
  IF n_elements(zsource) EQ 0 THEN zsource = .4
  IF n_elements(nlenses) EQ 0 THEN nlenses = 1
  IF NOT keyword_set(oplot) THEN oplot=0
  IF n_elements(sigrange) EQ 0 THEN sigrange = [600., 1400.]
  IF n_elements(SN0) EQ 0 THEN SN0  = 1. ;maximum is for SN0 = 1

  sigma = findgen(200)
  sigma = arrscl( sigma, sigrange[0], sigrange[1] )
  
  shapenoise=.5
  Dls = angdist(zsource, zlens)
  Ds  = angdist(zsource)

  x0 = 10.*(sigma/600.)^2*Dls/Ds ; Einstein radius
  
  n0 = 7500.                    ; sgal per square degree
  n0 = n0/3600.^2               ; per square arcsec (~2 square arcminute)
  n0 = n0*nlenses

  N = !pi*n0*x0^2               ;Number of galaxies in Einstein Radius

  S2N = 2.*sqrt(!pi)*sqrt(N)/shapenoise*SN0

  ytitle = 'S/N'
  xtitle = 'sigma'
  title = 'Isothermal Sphere.'

  IF n_elements(yrange) EQ 0 THEN yrange=[min(S2N),max(S2N)]
  IF NOT oplot THEN BEGIN
      plot, sigma, S2N, ytitle=ytitle, xtitle=xtitle, title=title,yrange=yrange
  ENDIF ELSE oplot, sigma, S2N


  return
END 
