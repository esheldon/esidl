PRO newphotoz_formula, meangauss, siggauss, modmean, modsig

  ;; Try new formula for redshift that is gaussian
  ;; for mean >> sigma but goes to zero as z^2 at low
  ;; z:  try modification factor:  z^2/( (sig/meanz)^2 + z^2 )

  nz = 1000
  minz = meangauss - 3.5*siggauss < 0.0
  maxz = meangauss + 4.5*siggauss
  z = arrscl( findgen(nz), minz, maxz )

  gaussz = gaussprob(z, meangauss, siggauss)

  w=where(z GT 0)
  modgauss = gaussz[w]*z[w]^2/( (siggauss/meangauss)^2 + z[w]^2)

  npts = 100
  norm = qgauss(modgauss, z[w], npts)
  modgauss = modgauss/norm

  ;; mean/var of new gauss
  modmean = qgauss(modgauss*z[w], z[w], npts)
  modvar  = qgauss(modgauss*(z[w]-modmean)^2, z[w], npts)
  modsig = sqrt(modvar)

  ;; calculate variance of new function
  print, meangauss, modmean
  print, siggauss,   modsig

  maxplot = 100
  xrange = [minz,maxz]
  plot, z[w], modgauss, xrange=xrange
  oplot, z[w], modgauss, color=!green
  oplot, z, gaussz
  oplot, [0, 0], [0, maxplot]

  oplot, [modmean, modmean], [0, maxplot], color=!red
  oplot, [modmean+modsig, modmean+modsig], [0, maxplot], color=!red, line=2
  oplot, [modmean-modsig, modmean-modsig], [0, maxplot], color=!red, line=2


  
END 
