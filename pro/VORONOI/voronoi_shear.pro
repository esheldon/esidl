PRO voronoi_shear, lcat, clr

  IF n_params() LT 2 THEN BEGIN
      print,'-Syntax: voronoi_shear, lcat, clr'
      return
  ENDIF 

  print,'Calculating the Voronoi density for each foreground galaxy'
  voronoi_density, lcat.ra, lcat.dec, dens

  run1=752
  run2=756
  rmin=10.
  rmax=600.
  binsize=39.

  outdir='/sdss4/data1/esheldon/GAL_GAL/DENSITY/'
  indir = '/sdss3/data1/corrected/corr752/1/combined/'

  nrand = 30000L
  n=n_elements(lcat)
  s=sort(dens)

  lowcat = lcat[s[0: n/2.]]
  low_s = sort(lowcat.ra)

  denscat=lcat[s[n/2.+1: n-1]]
  dens_s = sort(denscat.ra)
  
  addstr = 'dense_'
  objshear_fix, 752, 756, clr, rmin, rmax, binsize, outdir=outdir, lcat=denscat[dens_s], addstr=addstr,nrand=nrand,indir=indir
  
  addstr = 'low_'
  objshear_fix, 752, 756, clr, rmin, rmax, binsize, outdir=outdir, lcat=lowcat[low_s], addstr=addstr,nrand=nrand,indir=indir

  return
END 
