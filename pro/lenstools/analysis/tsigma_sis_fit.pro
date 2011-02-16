PRO tsigma_sis_fit, shstruct, sismass, sismasserr, sissigma2, sissigmaerr2, wuse=wuse

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: tsigma_sis_fit, shstruct, sismass, sismasserr, sissigma, sissigmaerr, wuse=wuse'
      return
  ENDIF 

  ;; fit mean shear in annuli to SIS

  IF n_elements(wuse) EQ 0 THEN BEGIN 
      nrad = n_elements(shstruct.rmax_act)
      wuse=lindgen(nrad)
  ENDIF ELSE BEGIN 
      nrad = n_elements(wuse)
  ENDELSE 

  sismass = fltarr(nrad)
  sismasserr = sismass
  sissigma2 = sismass
  sissigmaerr2 = sismass

  w=where( shstruct.rmax_act[wuse] GT 0., nw)

  IF nw NE 0 THEN BEGIN 

      w=wuse[w]

      rmin_factor = (1. + shstruct.rmin[w]/shstruct.rmax_act[w])
      ;; sigma is in Msolar/pc^2, convert radius to pc
      mfac = rmin_factor*!pi*(shstruct.rmax_act[w]*1000.)^2
      sismass[w] = mfac*shstruct.tsigma[w]
      sismasserr[w] = mfac*shstruct.tsigmaerr[w]

      ;; convert rmax to Mpc
      rmaxMpc = shstruct.rmax_act[w]/1000.
      sissigma2[w] = sismass[w]/7.31e8/rmaxMpc ;gives sigma_v^2 in km/s^2
      sissigmaerr2[w] = sismasserr[w]/7.31e8/rmaxMpc 

  ENDIF ELSE BEGIN
      print,'No good values of rmax_act > 0!'
      return
  ENDELSE 
  return
END 
