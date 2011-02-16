PRO colprint2,v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, $
              v15,v16,v17,v18,lun=lun,format=format,wuse_in=wuse_in

  ;;THIS VERSION PUTS - IN FOR MISSING ELEMENTS

  npar = n_params()
  IF npar LT 1 THEN BEGIN 
      print,'-Syntax: colprint,v1, v2, v3, ...,v18,lun=lun,format=format,wuse=wuse'
      return
  ENDIF 

  flag=0
  nbin = n_elements(v1)
  maxb = nbin
  FOR i=2, npar DO BEGIN
      nbinold = nbin
      istr = ntostr(i)
      com = 'nbin = max([nbin, n_elements(v'+istr+')])'
      tt=execute(com)
      com='maxb = [maxb, n_elements(v'+istr+')]'
      tt=execute(com)
      IF nbin NE nbinold THEN flag=1
  ENDFOR 
  
  ;; did the user input wuse?
  IF n_elements(wuse_in) EQ 0 THEN wuse=lindgen(nbin) ELSE BEGIN 
      wuse = wuse_in
      max_wuse = max(wuse)

      IF (max_wuse GT nbin) THEN BEGIN 
          message,'subscripts in wuse are too big: ignoring', /Informational
          wuse = lindgen(nbin)
      ENDIF ELSE nbin = n_elements(wuse)
  ENDELSE 

  IF flag THEN print,'Arrays not equal length.'

  IF n_elements(lun) EQ 0 THEN BEGIN
      com = 'print, ' 
  ENDIF ELSE BEGIN
      com='printf, '+ntostr(lun)+', '
  ENDELSE 
  IF n_elements(format) NE 0 THEN BEGIN
      endstr = ", format='"+format+"'"
  ENDIF ELSE BEGIN
      endstr=''
  ENDELSE 

  FOR iibin = 0, nbin-1 DO BEGIN 

      ibin = wuse[iibin]
      
      command = com

      binstr = ntostr(ibin)
      FOR i=1L, npar DO BEGIN 
      
          istr = ntostr(i)
          IF ibin GE maxb[i-1] THEN BEGIN 
              command = command + ' "         -   " '
          ENDIF ELSE BEGIN 
              command = command+'v'+istr+'['+binstr+']'
          ENDELSE 
          IF i NE npar THEN command = command+', '

      ENDFOR 
  
      command = command + endstr
      
      test=execute(command)
      IF test NE 1 THEN return

  ENDFOR 

return 
END 
  
