PRO voronoi_simpwtheta, lcat, clr

  IF n_params() LT 2 THEN BEGIN
      print,'-Syntax: voronoi_simpwtheta, lcat, clr'
      return
  ENDIF 

  print,'Calculating the Voronoi density for each foreground galaxy'
  voronoi_density, lcat.ra, lcat.dec, dens

  run1=752
  run2=756
  rmin=10.
  rmax=1200.                    ;Can go far out on wtheta
  binsize=50.

  nrand = 150000L

  outdir='/sdss4/data1/esheldon/TMP/'

  n=n_elements(lcat)
  s=sort(dens)

  lowcat = lcat[s[0: n/2.]]
  low_s = sort(lowcat.ra)

  denscat=lcat[s[n/2.+1: n-1]]
  dens_s = sort(denscat.ra)

  ;; use foreground sample as sources to get out excess galaxies.

  objshear, 752, 756, clr, rmin, rmax, binsize, outdir=outdir, lcat=denscat[dens_s], scat=denscat[dens_s], nrand=nrand
  objshear, 752, 756, clr, rmin, rmax, binsize, outdir=outdir, lcat=lowcat[low_s], scat=lowcat[low_s], nrand=nrand

  return
END 
