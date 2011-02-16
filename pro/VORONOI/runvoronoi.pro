PRO runvoronoi, lcat, stripe

  ;; run voronoi on tsgals_spec file, add voronoi_density tag
  ;; and overwrite file

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: runvoronoi, lcat, stripe'
      return
  ENDIF 

  ;; rotate ra,dec to the equator where it is 
  ;; rectangular

  stripe_inclination, stripe, inc
  ra=lcat.ra
  dec=lcat.dec
  rotate2equator, ra, dec, inc, ra_prime, dec_prime, /noresetdomain

  voronoi_density, ra_prime, dec_prime, voronoi_dens
  lcat.voronoi_dens = voronoi_dens

return
END 
