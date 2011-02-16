PRO plot_lens_density_contrast, dops=dops, encap=encap

  dir = '~/lensout/combstripe/comb/'
  file = dir + 'zgal_gal_stripe09_10_11_12_13_14_15_27_28_29_30_31_32_33_34_35_36_37_gri_recorr_h_jack_comb_comoving_N1.fit'

  struct = mrdfits(file, 1)

  psfile = '~/plots/deltasig'
  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN BEGIN 
          psfile = psfile+'.eps' 
          begplot,name=psfile, /encap
      ENDIF ELSE BEGIN 
          psfile=psfile+'.ps'
          begplot,name=psfile
      ENDELSE 
  ENDIF 

  xrange = [0.015, 17.0]
  yrange = [0.19, 130]
  plot_density_contrast,struct, /log, /axis_shear, /aspect, /delta, /mpc,$
                        xrange=xrange, yrange=yrange, /sigma_comoving
  oplot, struct.meanr/1000, struct.yfit*struct.meanscritinv*(1+struct.zmean)^2

  IF keyword_set(dops) THEN BEGIN 
      IF keyword_set(encap) THEN endplot, /trim_bbox $
      ELSE endplot
  ENDIF 

END 
