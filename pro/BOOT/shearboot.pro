PRO shearboot, lenses, nresample, shear, shearerr, ortho, orthoerr, $
               tshear, tshearerr, tortho, torthoerr

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: shearboot, lenses, nresample, shear, shearerr, ortho, orthoerr, tshear, tshearerr, tortho, torthoerr'
      print,' lenses is the lensum output from zobjshear, objshear'
      return
  ENDIF 

  nlenses = n_elements(lenses)
  nbin = n_elements(lenses[0].rsum)
  
  ;; don't use last bin!!!!
  ;; nbin = nbin-1

  shearsamp = fltarr(nresample, nbin)
  orthosamp = shearsamp
  tshearsamp = shearsamp
  torthosamp = shearsamp

  seed = long(systime(1))
  FOR samp=0L, nresample-1 DO BEGIN 

      sind = round( (nlenses-1)*randomu(seed, nlenses) )

      wsum = 0.
      tsh = 0.
      tor = 0.

      FOR bin=0L, nbin-1 DO BEGIN 
          wsum1 = total( lenses[sind].wsum[bin] )
          tsh1 =  total( lenses[sind].etansum[bin] )
          tor1 =  total( lenses[sind].eradsum[bin] )

          shearsamp[samp, bin] = tsh1/wsum1/2.
          orthosamp[samp, bin] = tor1/wsum1/2.

          wsum = wsum + wsum1
          tsh  = tsh  + tsh1
          tor  = tor  + tor1
          
          tshearsamp[samp, bin] = tsh/wsum/2.
          torthosamp[samp, bin] = tor/wsum/2.

      ENDFOR 
      IF samp MOD 100 EQ 0 THEN print,'.',format='(a,$)'
  ENDFOR 

  bootstrap, shearsamp, shear, shearerr
  bootstrap, orthosamp, ortho, orthoerr
  bootstrap, tshearsamp, tshear, tshearerr
  bootstrap, torthosamp, tortho, torthoerr

  orthosamp=0
  shearsamp=0
  torthosamp=0
  tshearsamp=0

return
END 

