PRO compare_lumavg_regavg, lave, dops=dops

  IF keyword_set(dops) THEN BEGIN 
      psfile_front = '~/plots/compare_lumavg_regavg'
  ENDIF 

  dir = esheldon_config("lensout_dir")+'combstripe/comb/'

  t=mrdfits(dir + 'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_N1.fit',1)

  ldir = dir + 'sublum/r/'

  l1=mrdfits(ldir+'lum1threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_N1.fit',1)
  l2=mrdfits(ldir+'lum2threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_N1.fit',1)
  l3=mrdfits(ldir+'lum3threebinnum_zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_N1.fit',1)

  lave = t
  nbin = n_elements(t.meanr)
  FOR i=0L, nbin-1 DO BEGIN 

      sig = [l1.sigma[i], l2.sigma[i], l3.sigma[i]]
      err = [l1.sigmaerr[i], l2.sigmaerr[i], l3.sigmaerr[i]]

      wmom, sig, err, wmean, wsig, werr

      lave.sigma[i] = wmean
      lave.sigmaerr[i] = werr

  ENDFOR 

  calc_r0_gamma_wgm,lave,$
                    /dolegend,yfit=yfit,wuse=wuse, $
                    /replace, $
                    noprompt=noprompt, r0range=[3.5,10.5],gamrange=[1.5,1.95],$
                    psfile_front=psfile_front,nr0=100,ngam=100

  

  IF keyword_set(dops) THEN begplot,name=psfile_front+'.ps',/color,/encap
  IF !d.name EQ 'X' THEN key=get_kbrd(1)
  IF !d.name EQ 'PS' THEN color=!blue ELSE color=!green
  plot_density_contrast,t,/mpc,/log,aspect=1,/delta
  oploterror, lave.meanr/1000., lave.sigma, lave.sigmaerr,psym=4,$
              color=color,errc=color
  
  oplot, t.meanr/1000,t.yfit,color=color
  oplot, lave.meanr/1000,lave.yfit

  legend,['Lum averaged', 'Standard'], psym=[8,4], color=[!p.color,color],$
         charsize=1,/right,box=0,thick=[!p.thick,!p.thick]

  IF keyword_set(dops) THEN endplot
END 
