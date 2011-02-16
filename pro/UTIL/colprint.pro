PRO colprint,v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, $
             v15,v16,v17,v18,lun=lun,format=format,file=file,wuse_in=wuse_in, $
             header=header

  npar = n_params()
  IF npar LT 1 THEN BEGIN 
      print,'-Syntax: colprint,v1, v2, v3, ...,v18,lun=lun,file=file,format=format, wuse=wuse'
      return
  ENDIF 

  dofile = 0
  IF n_elements(lun) EQ 0 THEN BEGIN 
      if n_elements(file) ne 0 then begin 
          dofile = 1
          openw, lun, file, /get_lun
      endif else begin
          lun = -1
      endelse
  ENDIF 

  flag=0
  nbin = n_elements(v1)
  FOR i=2L, npar DO BEGIN
      nbinold = nbin
      istr = ntostr(i)
      com = 'nbin = n_elements(v'+istr+')'
      tt=execute(com)
      IF nbin NE nbinold THEN flag=1
      nbin = min([nbin,nbinold])
  ENDFOR 
  
  IF flag THEN BEGIN
      print
      print,'*** Arrays not equal length. Not printing'
      print
      return
  ENDIF 

  ;; did the user input wuse?
  IF n_elements(wuse_in) EQ 0 THEN wuse=lindgen(nbin) ELSE BEGIN 
      wuse = wuse_in
      max_wuse = max(wuse)

      IF (max_wuse GT nbin) THEN BEGIN 
          message,'subscripts in wuse are too big: ignoring', /Informational
          wuse = lindgen(nbin)
      ENDIF ELSE nbin = n_elements(wuse)
  ENDELSE 

  comfront = 'printf, lun, '
  IF n_elements(format) NE 0 THEN BEGIN
      endstr = ", format='"+format+"'"
  ENDIF ELSE BEGIN
      endstr=''
  ENDELSE 

  if n_elements(header) ne 0 then printf, lun, header
  FOR iibin = 0L, nbin-1 DO BEGIN 

      ibin = wuse[iibin]

      command = comfront

      binstr = ntostr(ibin)
      FOR i=1L, npar DO BEGIN 
      
          istr = ntostr(i)
          command = command+'v'+istr+'['+binstr+'], '
      
      ENDFOR 
  
      command = command + 'format=format'
      
      test=execute(command)
      IF test NE 1 THEN return

  ENDFOR 

  IF dofile THEN BEGIN
      flush, lun
      free_lun, lun
  ENDIF 

return 
END 
  
