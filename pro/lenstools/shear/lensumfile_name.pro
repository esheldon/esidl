PRO lensumfile_name, stripes, front, clr, dir, name, $
                     hirata=hirata, recorr=recorr, $
                     lrg_sources=lrg_sources, rlrg_sources=rlrg_sources, $
                     subLumClr=subLumClr, $
                     meanFile=meanFile, $
                     ext=ext

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: lensumfile_name, stripes, front, clr, dir, name, $'
      print,'           /hirata, /recorr, $'
      print,'           /lrg_sources, /rlrg_sources, $'
      print,'           subLumClr=subLumClr, $'
      print,'           /meanFile, $'
      print,'           ext=ext'
      return
  ENDIF 

  nLumClr = n_elements(subLumClr)
  IF nLumCLR NE 0 THEN BEGIN 
      addDir = 'sublum/'+!colors[subLumClr]+'/'
  ENDIF ELSE BEGIN 
      addDir = ''
  ENDELSE  

  IF keyword_set(recorr) THEN recorrstr = '_recorr' ELSE recorrstr=''
  IF keyword_set(hirata) THEN hiratastr = '_h' ELSE hiratastr=''

  IF keyword_set(lrg_sources) THEN BEGIN 
      lrgstr = '_lrg' 
  ENDIF ELSE BEGIN 
      IF keyword_set(rlrg_sources) THEN lrgstr = '_rlrg' ELSE lrgstr=''
  ENDELSE 

  IF n_elements(ext) EQ 0 THEN ext='N1.fit'

  stripe_string = stripearr2string(stripes)
  clr_string = clrarr2string(clr)

  IF keyword_set(meanFile) THEN fTypeStr = '' ELSE fTypeStr='_lensum'
  filebase = $
    'stripe'+stripe_string+'_'+clr_string+lrgstr+recorrstr+hiratastr+$
    fTypeStr+'_'+ext

  dir = esheldon_config("lensout_dir") + 'stripe'+stripe_string+'/'+addDir
  name = front + filebase

END 
