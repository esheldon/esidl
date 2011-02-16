PRO spectra_combine_rand_lensout

  stripes = [9,10,11,12,13,14,15,$
             27,28,29,30,31,32,33,34,35,36,37,$
             76,86]
  randNum = 20 + lindgen(10)

  ;; All photoz
  combine_rand_lensout, randNum, type, stripes=stripes, /overwrite
  combine_rand_lensout, randNum, 'matchLRGNew_', stripes=stripes, /overwrite


  FOR subLumClr=2,2 DO BEGIN 

      FOR ilum=1,3 DO BEGIN 
          ltype = 'lum'+ntostr(iLum)+'threebinnum_matchLRGNew_'

          combine_rand_lensout, randNum, ltype, subLumClr=subLumClr, $
            stripes=stripes, /overwrite

      ENDFOR 

  ENDFOR 

  ;; rachel lrg photoz
  combine_rand_lensout, randNum, type, stripes=stripes, /overwrite, $
    /rlrg_sources, ext="N2.fit"
  combine_rand_lensout, randNum, 'matchAll_', stripes=stripes, /overwrite, $
    /rlrg_sources, ext="N2.fit"
  
  FOR subLumClr=2,2 DO BEGIN 

      FOR ilum=1,3 DO BEGIN 
          ltype = 'lum'+ntostr(iLum)+'threebinnum_matchAll_'

          combine_rand_lensout, randNum, ltype, subLumClr=subLumClr, $
            stripes=stripes, /overwrite, /rlrg_sources, ext="N2.fit"

      ENDFOR 

  ENDFOR 

return


  FOR subLumClr=0,4 DO BEGIN 

      FOR ilum=1,4 DO BEGIN 
          ltype = 'lum'+ntostr(iLum)+'fourbinnum_'

          combine_rand_lensout, randNum, ltype, subLumClr=subLumClr, $
            /overwrite

      ENDFOR 

  ENDFOR 



END 
