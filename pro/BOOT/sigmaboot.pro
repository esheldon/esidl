PRO sigmaboot, lenses, nresample, sigma, sigmaerr, orthosig, orthosigerr, $
               tsigma, tsigmaerr, torthosig, torthosigerr

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: sigmaboot, lenses, nresample, sigma, sigmaerr, orthosig, orthosigerr, tsigma, tsigmaerr, torthosig, torthosigerr'
      print,' lenses is the lensum output from zobjshear, objshear'
      return
  ENDIF 

  nlenses = n_elements(lenses)
  nbin = n_elements(lenses[0].rsum)
  
  ;; don't use last bin!!!!
  ;; nbin = nbin-1

  sigmasamp = fltarr(nresample, nbin)
  orthosigsamp = sigmasamp
  tsigmasamp = sigmasamp
  torthosigsamp = sigmasamp

  seed = long(systime(1))
  FOR samp=0L, nresample-1 DO BEGIN 

      sind = round( (nlenses-1)*randomu(seed, nlenses) )

      wsum = 0.
      tsi  = 0.
      tor  = 0.

      FOR bin=0L, nbin-1 DO BEGIN 
          wsum1 = total(lenses[sind].wsum[bin])
          tsi1  = total( lenses[sind].tansigsum[bin] )
          tor1  = total( lenses[sind].radsigsum[bin] )

          sigmasamp[samp, bin] = tsi1/wsum1/2.
          orthosigsamp[samp, bin] = tor1/wsum1/2.
          
          wsum = wsum + wsum1
          tsi  = tsi  + tsi1
          tor  = tor  + tor1

          tsigmasamp[samp, bin] = tsi/wsum/2.
          torthosigsamp[samp, bin] = tor/wsum/2.

      ENDFOR 
      IF samp MOD 100 EQ 0 THEN print,'.',format='(a,$)'
  ENDFOR 

  bootstrap, sigmasamp, sigma, sigmaerr
  bootstrap, orthosigsamp, orthosig, orthosigerr
  bootstrap, tsigmasamp, tsigma, tsigmaerr
  bootstrap, torthosigsamp, torthosig, torthosigerr

  orthosigsamp=0
  sigmasamp=0
  torthosigsamp=0
  tsigmasamp=0

return
END 

