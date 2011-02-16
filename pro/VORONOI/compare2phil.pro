PRO compare2phil, sigma, cutoff, zl, zs

  IF n_params() LT 2 THEN BEGIN
      print,'-Syntax: compare2phil, sigma, cutoff, zl, zs'
      return
  ENDIF 


  IF n_elements(zl) EQ 0 THEN zl=.15
  IF n_elements(zs) EQ 0 THEN zs=.4

  zs = .4
  zl = .15
  philzl = .172

  phil_sigcrit = 1./.392
  my_sigcrit = sigmacrit(zs, zl, h=1.)

  philDl = angdist(philzl, h=1.)
  myDl = angdist(zl, h=1.)
  

  rcut = ntostr(long(cutoff))
  indir = '/sdss4/data1/esheldon/TMP/conv/'

  myfile=indir+'conv_tot_r_cut'+rcut+'.txt'

  IF NOT exist(myfile) THEN BEGIN
      print,'No file ',myfile
      return
  ENDIF 

  readcol,myfile,radius,kappa,kappa_int,/silent
  myshear = (kappa_int - kappa)*(sigma)^2
  ;; rescale to phil's sigmacrit
  myshear = myshear*my_sigcrit/phil_sigcrit*myDl/philDl
  
  philfile = '/sdss3/usrdevel/esheldon/idl.lib/VORONOI/galphilshear.dat'
 
  readcol,philfile,meanr,shear,/silent

  plot,radius,myshear
  oplot,meanr,shear,psym=2

  outf = interpol(myshear,radius,meanr)
  oplot,meanr,outf,psym=1
  print,'     meanr     mine     phil'
  forprint,meanr,outf,shear
  

return
END
