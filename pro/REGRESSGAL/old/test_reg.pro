pro test_reg, lcat,scat

  IF (n_elements(lcat) EQ 0) OR (n_elements(scat) EQ 0) THEN BEGIN 
      
      infile='/sdss3/usrdevel/davej/sim-sg-data2.fits'
      t=mrdfits(infile,1)

      nb=n_elements(t.xb)
      nf=n_elements(t.xf)

      tlcat=create_struct('lambda', 0d, $
                         'eta',0d,$
                         'z1d', 0.0)

      tscat=create_struct('lambda', 0d, $
                         'eta', 0d, $
                         'momerr', 0.0, $
                         'e1', 0.0, $
                         'e2', 0.0)

      lcat = replicate(tlcat, nf)
      scat = replicate(tscat, nb)

      lcat.lambda = t.xf
      lcat.eta = t.yf
      lcat.z1d = t.zf


      
      scat.lambda = t.xb
      scat.eta = t.yb
      scat.e1 = t.e1
      scat.e2 = t.e2
      scat.momerr = 0.0

      ss=sort(scat.lambda)
      scat=scat[ss]

  ENDIF 

  rmax=3000.
  rmin=20.
  binsize=80.

  stripe=10
  clr=2

  step=100

  zregressgal, stripe, lcat, scat, clr, rmin, rmax, binsize, step=step,$
    /simulation,addstr='sim'


return
end


