PRO shapelet_st2fits, pattern=pattern

  dir = '~/shapelet_outputs/'
  IF n_elements(pattern) EQ 0 THEN BEGIN 
      files = file_search(dir+'stripe[0-9][0-9]-vagc-shapelets-*.st', count=count)
  ENDIF ELSE BEGIN
      files = file_search(dir+pattern, count=count)
  ENDELSE 

  FOR i=0L, count-1 DO BEGIN 
      print,files[i]
  ENDFOR 

END 
