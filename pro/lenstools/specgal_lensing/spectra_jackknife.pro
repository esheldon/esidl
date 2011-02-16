PRO spectra_jackknife

  randNum = lindgen(10)

  shear_jackknife
  shear_jackknife, randNum=randNum

  FOR subLumClr=0,4 DO BEGIN 

      FOR ilum=1,4 DO BEGIN 
          type = 'lum'+ntostr(iLum)+'fourbinnum_'

          shear_jackknife, type, subLumClr=subLumClr
          shear_jackknife, type, subLumClr=subLumClr, randNum=randNum
      ENDFOR 

  ENDFOR 


END 
