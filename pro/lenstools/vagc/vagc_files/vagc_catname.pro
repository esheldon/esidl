FUNCTION vagc_catname, lss=lss, letter=letter, post=post, sample=sample, $
                       dir=dir

  ;; Define the lss sample
  IF n_elements(sample) EQ 0 THEN sample = sdssidl_config('lss_vers')
  IF n_elements(letter) EQ 0 THEN letter = 'full'
  IF n_elements(post)   EQ 0 THEN post = '0'

  IF keyword_set(lss) THEN BEGIN 
      catname = 'lss_'+sample+letter+post
  ENDIF ELSE BEGIN 
      catname = 'vagc'
  ENDELSE 

  IF arg_present(dir) THEN dir=vagc_lensinput_dir(lss=lss)
  return, catname

END 
