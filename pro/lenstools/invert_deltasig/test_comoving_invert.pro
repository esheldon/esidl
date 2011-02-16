PRO test_comoving_invert, dops=dops
  
  ;; This is to make sure that the order of the comoving thing
  ;; doesn't matter.  First do all at z=zmean, then convert to
  ;; comoving and do again.  Then compare by converting all to
  ;; comoving

  IF n_elements(nrealize) EQ 0 THEN nrealize = 1

  name='~/plots/test_comoving_invert.ps'
  print,name
  
  ;; get radii from iro's points
  dir = '~/lensout/combstripe/comb/'
  dfile = dir +  'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_N1.fit'

  t = mrdfits(dfile, 1)

  zmean = t.zmean
  meanr = t.meanr/1000
  deltasig = t.sigma
  deltasig_err = t.sigmaerr
  covariance = t.covariance

  meanr_comoving = meanr*(1+zmean)
  deltasig_comoving = deltasig/(1+zmean)^2
  deltasig_err_comoving = deltasig_err/(1+zmean)^2
  covariance_comoving = covariance/(1+zmean)^4

  
  ;; dave's code
  print,'Doing redshift = '+ntostr(zmean)
  print,'Calling delta_sigma2der_sigma'
          
  delta_sigma2der_sigma, meanr, deltasig, dersig, covariance, $
                         covder, derr
  print,'Calling der_sigma2xi'
  rrange = [min(meanr), max(meanr)]
  der_sigma2xi, meanr, dersig, r3, xi, corrfac, covder, covxi, $
                rrange=rrange, zmean=zmean


  print,'Doing Comoving'
  print,'Calling delta_sigma2der_sigma'
          
  delta_sigma2der_sigma, meanr_comoving, deltasig_comoving, $
                         dersig_comoving, covariance_comoving, $
                         covder_comoving, derr_comoving
  print,'Calling der_sigma2xi'
  rrange_comoving = [min(meanr_comoving), max(meanr_comoving)]
  der_sigma2xi, meanr_comoving, dersig_comoving, $
                r3_comoving, xi_comoving, corrfac_comoving, $
                covder_comoving, covxi_comoving, $
                rrange=rrange_comoving

  !p.multi=0

  IF keyword_set(dops) THEN begplot,name=name,/color
 
  setup_mystuff

  ;; plots
  IF !d.name EQ 'PS' THEN oclr = !blue ELSE oclr = !green

  aplot, 1, r3*(1+zmean), xi, psym=8, /xlog, /ylog
  oplot, r3_comoving, xi_comoving, psym=4, color=oclr

  legend, ['Method 1','Method 2'], psym=[8,4], color=[!p.color,oclr],/right,$
          box=0, thick=[!p.thick,!p.thick]

  IF keyword_set(dops) THEN endplot

return
END 
