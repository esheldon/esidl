PRO srcgal_pixelmaskfile_name, stripe, clr, file, hirata=hirata, rlrgMask=rlrgMask

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: srcgal_pixelmaskfile_name, stripe, clr, file, hirata=hirata, /rlrgMask'
      return
  ENDIF 

  maskDir = sdssidl_config('shapecorr_dir') + 'masks/'
  clr_string = clrarr2string(clr)
  stripe_string = stripearr2string(stripe)

  IF keyword_set(hirata) THEN hirstr='_h' ELSE hirstr=''
  IF keyword_set(rlrgMask) THEN rlrgStr='_rlrg' ELSE rlrgStr = ''

  file = $
    maskDir + 'stripe'+stripe_string+'_srcgal_'+clr_string+hirstr+$
    rlrgStr+'_pixel_maskflags.fit'

END 
