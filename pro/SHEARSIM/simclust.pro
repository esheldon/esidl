PRO simclust, inclust, sigma, incat=incat, outcat=outcat, shifty=shifty, $
              shiftx=shiftx

  IF n_params() EQ 0 THEN BEGIN 
      print,'-Syntax: simclust, clust, sigma, incat=incat'
      return
  ENDIF 

  r1=752
  r2=756
  clr=2
  slength = 120.
  rfac = 10.

  clust = inclust
  IF n_elements(shiftx) NE 0 THEN BEGIN 
      clust.dec = clust.dec + shiftx/3600.
  ENDIF 
  IF n_elements(shifty) NE 0 THEN BEGIN 
      clust.ra = clust.ra + shifty/3600.
  ENDIF 
  

  IF n_elements(incat) EQ 0 THEN BEGIN 

      kappa_map, r1, r2, clr, clust, slength=slength, rfac=rfac, $
                 scat=scat, wsource=wsource
      incat = scat[wsource]

  ENDIF 

  zsource = .4
  zlens = clust.z
  sis_shear, zlens, zsource, sigma, 1, 1, outcat, incat=incat, $
             incenter = [clust.dec, clust.ra ], /error

  kappa_map, r1, r2, clr, clust, slength=slength, rfac=rfac, scat=outcat

  return
END 
