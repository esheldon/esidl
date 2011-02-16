PRO huan_goods_addnoise_scripts, scriptnum

  ;; for 0.7 use 0.675
  ;; for 0.8 use 0.78
  ;; for 0.9 use 0.87

  ;; This time I skipped 0.87 since I had already done it

  CASE scriptnum OF
      1: seeing = [ 0.20, 0.25, 0.30, 0.35, 0.40, 0.45 ]
      2: seeing = [ 0.50, 0.55, 0.60, 0.65, 0.70, 0.76 ]
      3: seeing = [ 0.80, 0.85, 0.90, 0.95, 1.00, 1.05 ]
      4: seeing = [ 0.78, 0.675,1.10, 1.15, 1.20 ]
      ELSE: message,'Unknown script number: '+ntostr(scriptnum)
  ENDCASE 

  scriptDir = '/net/cheops2/home/esheldon/idlscripts/goods_addnoise/'
  scriptName = scriptDir + 'goods_addnoise'+ntostr(scriptnum)+'.sh'

  nseeing = n_elements(seeing)

  openw, lun, scriptName, /get_lun

  printf,lun,'#!/bin/bash'
  printf,lun

  FOR i=0L, nseeing-1 DO BEGIN 

      sstr = ntostr(seeing[i])

      printf,lun,'idl<<EOF'
      printf,lun,'  huan_goods_addnoise, "des5yr", seeing='+sstr
      printf,lun,'EOF'
      printf,lun
  ENDFOR 

  free_lun, lun
  spawn,['chmod','755',scriptName],/noshell

END 
