PRO s2n_analytic, zlens, zsource

  density = 7500./(3600.^2)     ; per square arcsec
  err = .5

  Ds = angdist(zsource)
  Dls = angdist(zsource, zlens)

  sigma = [500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400.]

  x0 = 10.*(sigma/600.)^2*Dls/Ds

  fac = sqrt(!pi*density*x0^2)/.5

  schneider = 2.                ;For nu1 = .05

  mine = 2.                     ;For R=infinity

  plot, sigma, schneider*fac
  oplot, sigma, schneider*fac, psym=1
  oplot, sigma, mine*fac
  oplot, sigma, mine*fac, psym=2

  legend, ['Mine', 'Schneider'], psym=[2,1]

  return
END 
