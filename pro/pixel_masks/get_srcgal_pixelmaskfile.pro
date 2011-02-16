PRO get_srcgal_pixelmaskfile, stripe, clr, maskStruct, status=status, silent=silent, rlrgMask=rlrgMask, hirata=hirata

  status=1
  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: get_srcgal_pixelmaskfile, stripe, clr, maskStruct, status=status, /silent, /rlrgMask'
      return
  ENDIF 

  srcgal_pixelmaskfile_name, stripe, clr, file, rlrgMask=rlrgMask, $
    hirata=hirata
  IF NOT fexist(file) THEN BEGIN
      delvarx, maskStruct
      message,'File does not exist: '+file
  ENDIF 
  
  IF NOT keyword_set(silent) THEN BEGIN 
      print
      print,'Reading pixelMaskFlags file: ',file
  ENDIF 
  maskStruct = mrdfits(file,1, silent=silent)

  IF datatype(maskStruct) EQ 'STC' THEN status=0 ELSE delvarx, maskStruct

END 
