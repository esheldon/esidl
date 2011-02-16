PRO lumdensboot, lenses, random, Nper, nresample, lumdens, lumdenserr, covariance,  lumdiff_resamp, wuse=wuse

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: lumdensboot, lenses, random, Nper, nresample, lumdens, lumdenserr, covariance, lumdiff_resamp, wuse=wuse'
      return
  ENDIF 

  tt=systime(1)

  nlenses = n_elements(lenses)
  IF n_elements(wuse) NE 0 THEN BEGIN 
      nbin = n_elements(wuse)
  ENDIF ELSE BEGIN 
      nbin = n_elements(lenses[0].rsum)
      wuse = lindgen(nbin)
  ENDELSE   
  print
  print,'Number of objects: ',nlenses
  print,'Number of bins used: ',nbin
  print,'N resample: ',nresample
  lumdensamp = fltarr(nresample, nbin)
  rlumdensamp = lumdensamp

  sstr= ntostr(nresample)

  ;; Use histogram to get random for each lens
  ;; using histogram this way, with min=0, forces a bin for every
  ;; integer from 0 to the max. That way, we can use
  ;; arr1[i] as subscript for reverse_indices! 
  
  ;; match multi gets all matches, but also returnes reverse_indices
  ;; for our use.  Assumes each lens zindex is represented in random
  match_multi, lenses.zindex, random.zindex, rmatch, reverse_indices=revind

  srand = lonarr(Nper*Nlenses)
  ind = lonarr(Nper)

  print
  print,'Building Bootstrap samples'
  seed = long(systime(1))
  FOR samp=0L, nresample-1 DO BEGIN 

      sind = round( (nlenses-1)*randomu(seed, nlenses) )

      ;; get the randoms that match
      ;; This is faster than using match2rand again
      ;; by about a factor of 2
      ;; still small relative to the loop below, however
      ;; For Nbin=10, Nlens=16839, Nper=10 that loop takes
      ;; 20 times longer than this one
      ttt=systime(1)
      FOR iL=0L, Nlenses-1 DO BEGIN 
          zi = lenses[sind[iL]].zindex
          ;; use * to gaurantee we get Nper; it will crash otherwise
          ind[*] = revind[ revind[zi]:revind[zi+1] -1 ]
          srand[ iL*Nper: (iL+1)*Nper-1 ] = ind
      ENDFOR 
      print,systime(1)-ttt

      wsum = 0.
      ld = 0.
      rwsum = 0.
      rld = 0.

      ;; all the time is spent in this loop
      ttt=systime(1)
      FOR bin=0L, nbin-1 DO BEGIN 
          wbin = wuse[bin]

          wsum = total( lenses[sind].wsum[wbin] )
          ld =  total( lenses[sind].lsum[wbin] )

          rwsum = total( random[srand].wsum[wbin] )
          rld  = total( random[srand].lsum[wbin] )

          lumdensamp[samp, bin] = ld/wsum
          rlumdensamp[samp, bin] = rld/rwsum

      ENDFOR 
      print,systime(1)-ttt
      print,ntostr(samp+1)+'/'+sstr

      ;;IF (samp+1) MOD 10 EQ 0 THEN print,'.',format='(a,$)'
      ;;IF ( (samp+1) MOD 10 ) EQ 0 THEN print,ntostr(samp+1)+'/'+sstr

  ENDFOR 

  print
  print,'Bootstrapping'
  lumdiff_resamp = temporary(lumdensamp) - temporary(rlumdensamp)

  bootstrap, lumdiff_resamp, Nresample, Nbin, lumdens, lumdenserr, covariance

  ptime,systime(1)-tt

return
END 

