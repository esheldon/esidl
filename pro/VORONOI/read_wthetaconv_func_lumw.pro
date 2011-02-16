PRO read_wthetaconv_func_lumw, clr, wclr, beta, sigmav, radiuskpc, $
                               model_density_contrast, $
                               neighbor_density_contrast, $
                               central_density_contrast, $
                               central_density, $
                               model_ratio, model_ratio_cum, $
                               sigonly=sigonly, noweight=noweight


  IF n_params() EQ 0 THEN BEGIN 
      print,'-Syntax: read_wthetaconv_func, clr, wclr, beta, sigmav, '
      print,'         modelrkpc, model_density_contrast, '
      print,'         neighbor_density_contrast, central_density_contrast, '
      print,'         central_density, model_ratio, model_ratio_cum, '
      print,'         sigonly=sigonly, noweight=noweight'
      print,' type=0 (high dens) type=1 (low dens) type=2 (tot)'
      return
  ENDIF

  typestr = ['lum','lum1', 'lum2','lum3','lum4']

  setup_mystuff
  
  betastr = '0-2'
  indir = '/sdss5/data0/lensout/wtheta_conv/'
  type=0
  CASE 1 OF
      keyword_set(sigonly): savefile = indir + 'conv_'+!colors[wclr]+'w_'+$
        !colors[clr]+'_'+typestr[type]+'_sigonly_beta0-2.sav'
      keyword_set(noweight): savefile = indir + 'conv_'+!colors[wclr]+'w_'+$
        !colors[clr]+'_'+typestr[type]+'_noweight_beta0-2.sav'
      ELSE: savefile = indir + 'conv_'+!colors[wclr]+'w_'+!colors[clr]+'_'+typestr[type]+'_beta'+betastr+'.sav'
  ENDCASE 
  print,'Reading: ',savefile

  IF fexist(savefile) THEN restore, savefile $
  ELSE message,'File does not exist: '+savefile

END 
