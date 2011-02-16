PRO testzdist


  gfile = '/sdss3/usrdevel/esheldon/idl.lib/SIGMA_CRIT/coeffg.fit'
  rfile = '/sdss3/usrdevel/esheldon/idl.lib/SIGMA_CRIT/coeffr.fit'
  ifile = '/sdss3/usrdevel/esheldon/idl.lib/SIGMA_CRIT/coeffi.fit'


  coeffg = mrdfits(gfile)
  coeffr = mrdfits(rfile)
  coeffi = mrdfits(ifile)

  nn = 1123
  zmin = .01
  zmax = 1.5
  z = arrscl( findgen(nn), zmin, zmax)
  
  gzdist = fltarr(nn)
  rzdist = gzdist
  izdist = gzdist

  theta = .1

  FOR i=0, nn-1 DO BEGIN 

      f = (theta/z[i])^2 * exp(theta/z[i])/( exp(theta/z[i]) - 1 )^2

      arg = (z[i]-coeffg(1,*))^2 / (2*coeffg(2,*)^2) < 10.8
      gzdist[i] = total( coeffg(0,*) * exp(-arg )*f )

      arg = (z[i]-coeffr(1,*))^2 / (2*coeffr(2,*)^2)  < 10.8
      rzdist[i] = total( coeffr(0,*) * exp(-arg )*f )

      arg = (z[i]-coeffi(1,*))^2 / (2*coeffi(2,*)^2)  < 10.8
      izdist[i] = total( coeffi(0,*) * exp(-arg )*f )

  ENDFOR 

  gmean = total( gzdist*z )/total(gzdist)
  rmean = total( rzdist*z )/total(rzdist)
  imean = total( izdist*z )/total(izdist)

  plot, z, rzdist, line=0
  oplot, z, gzdist, line=1
  oplot, z, izdist, line=2
  yy = 1.e4

  legend, ['g-band <z> = '+ntostr(gmean), $
           'r-band <z> = '+ntostr(rmean), $
           'i-band <z> = '+ntostr(imean)], line=[1,0,2], /right

  key = get_kbrd(1)

  plot, z, rzdist/total(rzdist), line=0
  oplot, z, gzdist/total(gzdist), line=1
  oplot, z, izdist/total(izdist), line=2

  legend, ['g-band <z> = '+ntostr(gmean), $
           'r-band <z> = '+ntostr(rmean), $
           'i-band <z> = '+ntostr(imean)], line=[1,0,2], /right


  print,gmean,rmean,imean

return
END 
