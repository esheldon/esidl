PRO make_icons, expandstring

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: make_icons, expandstring'
      return
  ENDIF 

  ;; run in CWD

  files1 = findfile('*.jpg')
  files2 = findfile('*.JPG')
  files3 = findfile('*.gif')
  files4 = findfile('*.GIF')
  files5 = findfile('*.JPEG')
  files6 = findfile('*.jpeg')
  
  IF files1[0] NE '' THEN add_arrval, files1, imfiles
  IF files2[0] NE '' THEN add_arrval, files2, imfiles
  IF files3[0] NE '' THEN add_arrval, files3, imfiles
  IF files4[0] NE '' THEN add_arrval, files4, imfiles
  IF files5[0] NE '' THEN add_arrval, files5, imfiles
  IF files6[0] NE '' THEN add_arrval, files6, imfiles
  nfile = n_elements(imfiles)
  IF nfile EQ 0 THEN BEGIN
      print,'No files found'
      return
  ENDIF 

  FOR i=0L, nfile-1 DO BEGIN 
      
      xvcomm = 'xv -expand ' + expandstring+ ' '+imfiles[i]
      spawn, xvcomm

  ENDFOR 
  return
END 
