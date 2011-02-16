PRO netscape_image, imagedir

  ;; e.g. make_icons, dir, '.1'

  files1 = findfile(imagedir+'*.jpg')
  files2 = findfile(imagedir+'*.JPG')
  files3 = findfile(imagedir+'*.gif')
  files4 = findfile(imagedir+'*.GIF')
  files5 = findfile(imagedir+'*.JPEG')
  files6 = findfile(imagedir+'*.jpeg')
      

  IF files1[0] NE '' THEN add_arrval, files1, imfiles
  IF files2[0] NE '' THEN add_arrval, files2, imfiles
  IF files3[0] NE '' THEN add_arrval, files3, imfiles
  IF files4[0] NE '' THEN add_arrval, files4, imfiles
  IF files5[0] NE '' THEN add_arrval, files5, imfiles
  IF files6[0] NE '' THEN add_arrval, files6, imfiles
  IF n_elements(imfiles) EQ 0 THEN BEGIN
      print,'No files found'
      return
  ENDIF 
  
  nf=n_elements(imfiles)

  FOR i=0L, nf-1 DO BEGIN 
      
      command = "netscape -remote 'openFile("+imfiles[i]+")'"
      spawn,command,result
      IF (nf GT 1) AND (i NE nf-1) THEN BEGIN
          print,'Next object (y/n)?'
          key=get_kbrd(1)
          IF (key EQ 'n') OR (key EQ 'N') OR (key EQ 'q') THEN return
      ENDIF 
  ENDFOR 

return
END 
