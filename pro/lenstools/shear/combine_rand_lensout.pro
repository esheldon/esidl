PRO combine_rand_lensout, randNum, type, ext=ext, $
                          stripes=stripes, clr=clr, $
                          lrg_sources=lrg_sources, rlrg_sources=rlrg_sources, $
                          overwrite=overwrite, subLumClr=subLumClr, $
                          MaxBCG=MaxBCG

  IF n_elements(type) EQ 0 THEN type=''
  IF n_elements(ext) EQ 0 THEN ext='N1.fit'

  IF n_elements(stripes) EQ 0 THEN BEGIN 
      stripes = [9,10,11,12,13,14,15,$
                 27,28,29,30,31,32,33,34,35,36,37,$
                 76,86]
  ENDIF 
  IF n_elements(clr) EQ 0 THEN BEGIN 
      clr = [1,2,3]
  ENDIF 

  nRand = n_elements(randNum)

  IF keyword_set(MaxBCG) THEN BEGIN 
      front = type + 'MaxBCG_zrand'+ntostr(randNum)+'_'
      outFront = type + 'MaxBCG_zrand_'
  ENDIF ELSE BEGIN 
      front = type + 'zrand'+ntostr(randNum)+'_'
      outFront = type + 'zrand_'
  ENDELSE 

  nFile = n_elements(front)

  lensumfile_name, stripes, outFront, clr, dir, outFile, $
    lrg_sources=lrg_sources, rlrg_sources=rlrg_sources, $
    /meanFile, $
    /hirata, /recorr, ext=ext, subLumClr=subLumClr

  outFile = dir + outFile
  print
  print,'Will output to file: ',outFile
  print

  rFiles = strarr(nFile)
  FOR i=0L, nfile-1 DO BEGIN 
      lensumfile_name, stripes, front[i], clr, dir, tFile, $
        lrg_sources=lrg_sources, rlrg_sources=rlrg_sources, $
        /meanFile, $
        /hirata, /recorr, ext=ext, subLumClr=subLumClr

;      print,tFile
      rFiles[i] = dir + tFile
  ENDFOR 

  print
  print,'Combining files'
  combine_zshear, rFiles, shstruct

  shhdr = zshhdr(shstruct)
  IF NOT fexist(outFile) OR keyword_set(overwrite) THEN BEGIN  
      print
      print,'Writing to file: ',outFile
      mwrfits, shstruct, outFile, shhdr, /create
  ENDIF ELSE BEGIN 
      print
      print,'File already exists. Set /overwrite...'
      return
  ENDELSE 
END 
