FUNCTION vagc_lensinput_dir, lss=lss
  
  combdir = sdssidl_config('SHAPECORR_DIR') + 'combined/'
  IF keyword_set(lss) THEN BEGIN 
      dir = combdir + 'lss/'
  ENDIF ELSE BEGIN 
      dir = combdir + 'vagc/'
  ENDELSE 

  return,dir

END 
