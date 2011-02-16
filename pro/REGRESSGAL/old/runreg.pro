PRO runreg,clr,lcat,scat,FIRSTRA, LASTRA

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: runreg, clr [, lcat, scat,FIRSTRA, LASTRA'
      return
  ENDIF 

  run1 = 752
  run2 = 756
  r1str=ntostr(run1)
  r2str=ntostr(run2)

  colors=['u','g','r','i','z']

  binsize = 39.
  rmin = 10.
  rmax = 600.

  dir='/sdss4/data1/esheldon/CORRECTED/'
  lf=dir+'run'+r1str+'_'+r2str+'_lensgal_'+colors[clr]+'_overlap.fit'
  sf=dir+'run'+r1str+'_'+r2str+'_srcgal_'+colors[clr]+'_overlap.fit'
  
  lcat=mrdfits(lf,1)
  scat=mrdfits(sf,1)

  firstra=min(scat.ra)
  lastra=max(scat.ra)

  regressgal, run1, run2, clr, rmin, rmax, binsize, scat=scat, lcat=lcat, $
    FIRSTRA=FIRSTRA, LASTRA=LASTRA

  return
END 
