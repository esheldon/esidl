PRO mk_image_html_noicon, imagedir, htmlname, comments=comments, filelist=filelist

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: mk_image_html, imagedir, htmlname,comments=comments'
      return
  ENDIF 

  IF NOT keyword_set(comments) THEN comments=0

  htmlfile = imagedir+htmlname

  openw, lun, htmlfile, /get_lun

  printf, lun, '<HTML>'
  printf, lun, '<HEAD><TITLE>'+htmlname+'</TITLE></HEAD>'
  printf, lun, '<BODY link="#0066ff" vlink="#FF0000" bgcolor="#000000" text="#ffffff">'

  IF NOT keyword_set(filelist) THEN BEGIN 

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

  ENDIF ELSE BEGIN 
      imfiles = filelist
  ENDELSE 


  nfile = n_elements(imfiles)


  cmmt = ' '
  FOR i=0, nfile-1 DO BEGIN 

      IF comments THEN BEGIN 
          xvcomm = 'xv '+imfiles[i]+' &'
          spawn,xvcomm
          cmmt = ' '
          str = '"Comments for '+imfiles[i]+'"'
          print,format='($, '+str+')'
          read, cmmt

          printf,lun,cmmt+'<br>'
      ENDIF ELSE BEGIN 
          printf,lun,imfiles[i]+'<br>'
      ENDELSE 

      printf,lun,'<img src="./'+imfiles[i]+'"><br><br>'

  ENDFOR 

  printf,lun, '</BODY>'
  printf,lun, '</HTML>'

  free_lun, lun

return 
END 
