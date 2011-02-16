PRO rdlensout, struct, front=front, subLumClr=subLumClr, ext=ext, $
               silent=silent, meanFile=meanFile

  IF n_params() LT 1 THEN BEGIN  
      print,'-Syntax: rdlensout, struct, front=front, subLumClr=subLumClr, ext=ext, /silent, /meanFile'
      return
  ENDIF 

  IF n_elements(front) EQ 0 THEN front='zgal_gal_'

  dir = '~/lensout/'

  stripes = [9,10,11,12,13,14,15,$
             27,28,29,30,31,32,33,34,35,36,37,$
             76,86]
  clr = [1,2,3]

  lensumFile_name, stripes, front, clr, dir, name, ext=ext, $
    /hirata, /recorr, subLumClr=subLumClr, meanFile=meanFile

  file = dir + name

  IF NOT keyword_set(silent) THEN BEGIN 
      print
      print,'Reading file: ',file
  ENDIF 

  struct = mrdfits(file, 1, silent=silent)

return
  stripeString = stripearr2string(stripes)

  ss = 'stripe'+stripeString
 
  
  IF n_elements(ext) EQ 0 THEN ext='N1.fit'

  IF n_elements(front) NE 0 THEN BEGIN 
      file = dir + ss + '/'+front+'_'+ss+'_gri_recorr_h_lensum_'+ext
  ENDIF ELSE BEGIN 
      file = dir + ss + '/zgal_gal_'+ss+'_gri_recorr_h_lensum_'+ext
  ENDELSE 

  IF NOT keyword_set(silent) THEN BEGIN 
      print
      print,'Reading file: ',file
  ENDIF 

  struct = mrdfits(file, 1, silent=silent)

return
  mrdfits_multi, [f1, f2], struct, /diff,silent=silent

END 
