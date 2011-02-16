PRO read_wthetaconv_func, type, cutoff, sigmav, radiuskpc, model_density_contrast, $
                          neighbor_density_contrast, $
                          central_density_contrast, $
                          central_density, $
                          model_ratio, model_ratio_cum

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: read_wthetaconv_func, type, cutoff, sigmav, '
      print,'         modelrkpc, model_density_contrast, '
      print,'         neighbor_density_contrast, central_density_contrast, '
      print,'         central_density, model_ratio, model_ratio_cum'
      print,' type=0 (high dens) type=1 (low dens) type=2 (tot)'
      return
  ENDIF

  cutstr = '50-1200'
  typestr = ['high','low', 'tot']
  indir = '/sdss5/data0/lensout/wtheta_conv/'
  savefile = indir + 'conv_'+typestr[type]+'_cut'+cutstr+'.sav'


  IF fexist(savefile) THEN restore, savefile $
  ELSE message,'File does not exist: '+savefile

END 
