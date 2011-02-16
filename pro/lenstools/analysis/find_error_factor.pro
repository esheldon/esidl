PRO find_error_factor

  ;; calculate the factor that should be applied to error
  ;; bars on combined plots

  gg = 3.50959599
  rr = 2.235299689
  ii = 2.549197961
  gr = 1.586466533 
  gi = 1.6983151 
  ri = 1.512359722

  indir='/sdss5/data0/lensout/stripe10/'

  fcombcorr=indir+'main_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_comb_corr_N1.fit'
  fgcorr=indir+'main_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_g_corr_N1.fit'
  frcorr=indir+'main_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_r_corr_N1.fit'
  ficorr=indir+'main_zgal_gal_stripe10_stripe36_stripe37_stripe42_stripe43_stripe82_i_corr_N1.fit'

  combcorr=mrdfits(fcombcorr,1,/silent)
  gcorr=mrdfits(fgcorr,1,/silent)
  rcorr=mrdfits(frcorr,1,/silent)
  icorr=mrdfits(ficorr,1,/silent)

  scale = (1./3.)*( gcorr.shearerr^2/gg + rcorr.shearerr^2/rr + icorr.shearerr^2/ii )
  
  shear_bar2 = (gcorr.shearerr^2 + rcorr.shearerr^2 + icorr.shearerr^2 )/scale/3.
  error_fac = sqrt( (gr + gi + ri)/shear_bar2 )

  print,error_fac

;  erase & multiplot,[1,4]

;  yrange=[-10,40]

;  ploterror,gcorr.meanr,gcorr.sigma,gcorr.sigmaerr,psym=1,yrange=yrange
;  multiplot
;  ploterror,rcorr.meanr,rcorr.sigma,rcorr.sigmaerr,psym=1,yrange=yrange
;  multiplot
;  ploterror,icorr.meanr,icorr.sigma,icorr.sigmaerr,psym=1,yrange=yrange

;  multiplot
  ploterror,combcorr.meanr,combcorr.sigma,combcorr.sigmaerr*error_fac,psym=1,yrange=yrange
  oploterror,combcorr.meanr,combcorr.sigma,combcorr.sigmaerr*error_fac,psym=1,$
    color=!red,errcolor=!red
  oploterror,combcorr.meanr,combcorr.sigma,combcorr.sigmaerr,psym=1
;  multiplot,/reset


END 

