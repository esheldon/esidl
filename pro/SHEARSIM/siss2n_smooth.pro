PRO sis_s2n, sigma, zlens, zsource, n0=n0,oplot=oplot, yrange=yrange,nlenses=nlenses

;; finds s2n for the smoothing function I used

  IF n_elements(zsource) EQ 0 THEN zsource = .4
  IF n_elements(zlens)   EQ 0 THEN zlens   = .15
  IF n_elements(sigma)   EQ 0 THEN sigma   = 160 ;km/s
  IF n_elements(nlenses) EQ 0 THEN nlenses = 1
  IF NOT keyword_set(oplot) THEN oplot = 0

  shapenoise=.55
  Dls = angdist(zsource, zlens)
  Ds  = angdist(zsource)

  x0 = 10.*(sigma/600.)^2*Dls/Ds ; Einstein Radius

  n0 = 108000.
;  n0 = 7500.                    ; sgal per square degree
  n0 = n0/3600.^2               ; per square arcsec (~2 square arcminute)
  n0 = n0*nlenses

  N = !pi*n0*x0^2               ;Number of galaxies in Einstein Radius

  RS = findgen(200)             ; (R/s)
  RS = arrscl(RS, 6.5, 60.)

  SN0 = (1. - sqrt(2./!pi)/RS )/sqrt( 1. - 4./RS^2 )

  S2N = 2.*sqrt(!pi)*sqrt(N)/shapenoise*SN0;/sqrt(!pi)

  ytitle = 'S/N'
  xtitle = 'R/S'
  title = 'Isothermal Sphere  Sigma = '+ntostr(sigma)+$
           '  nlenses = '+ntostr(nlenses)

  max = max(S2N)
  IF n_elements(yrange) EQ 0 THEN yrange = [.9*max, 1.1*max]
  IF NOT oplot THEN BEGIN 
      
      ytitle= '(S/N)o'
      plot, RS, SN0, xtitle=xtitle, ytitle=ytitle, yrange=[.9,1.0]
  ENDIF 
  IF NOT oplot THEN BEGIN
      key=get_kbrd(1)
      plot, RS, S2N, ytitle=ytitle, xtitle=xtitle, title=title, yrange=yrange
  ENDIF ELSE oplot,RS,S2N


  return
END 
