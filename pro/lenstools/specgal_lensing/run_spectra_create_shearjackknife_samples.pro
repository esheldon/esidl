PRO run_spectra_create_shearjackknife_samples


;  spectra_create_shearjackknife_samples
;  spectra_create_shearjackknife_samples, randNum=lindgen(10)

  FOR subLumClr=0,4 DO BEGIN 
      FOR ilum=1,4 DO BEGIN 
          type = 'lum'+ntostr(iLum)+'fourbinnum_'

          spectra_create_shearjackknife_samples, type, subLumClr=subLumClr
          spectra_create_shearjackknife_samples, type, subLumClr=subLumClr,$
            randNum=lindgen(10)
      ENDFOR 
  ENDFOR 

END 
