PRO lumdensjackknife, lenses, random, lumdens, lumdenserr, covariance, wuse=wuse

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: lumdensjackknife, lenses, random, lumdens, lumdenserr, covariance, wuse=wuse'
      return
  ENDIF 

  Nlenses = n_elements(lenses)
  IF n_elements(wuse) NE 0 THEN BEGIN 
      Nbin = n_elements(wuse)
  ENDIF ELSE BEGIN 
      Nbin = n_elements(lenses[0].rsum)
      wuse = lindgen(nbin)
  ENDELSE  

  lumdens  = fltarr(Nbin)
  lumdenserr = lumdens
  covariance = fltarr(Nbin, Nbin)

  factor = (Nlenses - 1.)/Nlenses

  sampval = fltarr(Nlenses, Nbin)
  jval = fltarr(Nbin)
  Lstot = 0d
  Lsmod = 0d
  Lwtot = 0d
  Lwmod = 0d
  Lsmean = 0d


  Rstot = 0d
  Rsmod = 0d
  Rwtot = 0d
  Rwmod = 0d
  Rsmean = 0d

  smean = 0d

  print
  print,'Number of objects: ',Nlenses
  print,'Jackknifing'

  ;; Use histogram to get random for each lens
  ;; using histogram this way, with min=0, forces a bin for every
  ;; integer from 0 to the max. That way, we can use
  ;; arr1[i] as subscript for reverse_indices! 
  
  ;; match multi gets all matches, but also returnes reverse_indices
  ;; for our use.  Assumes each lens zindex is represented in random
  match_multi, lenses.zindex, random.zindex, rmatch, reverse_indices=revind
  
  ;; this gets all matches

  FOR bin=0L, Nbin-1 DO BEGIN 

      Lstot = total( lenses[*].lsum[bin], /double)
      Lwtot = total( lenses[*].wsum[bin], /double)
      Lsmean = Lstot/Lwtot
      
      Rstot = total( random[rmatch].lsum[bin], /double)
      Rwtot = total( random[rmatch].wsum[bin], /double)
      Rsmean = Rstot/Rwtot

      smean = Lsmean - Rsmean

      FOR j=0L, Nlenses-1 DO BEGIN 
          ;; lenses

          zi = lenses[j].zindex

          Lsmod = Lstot - lenses[zi].lsum[bin]
          Lwmod = Lwtot - lenses[zi].wsum[bin]
          
          ;; Random
          ind = revind[ revind[zi]:revind[zi+1] -1 ]
          Rsmod = Rstot - total( random[ind].lsum[bin], /double)
          Rwmod = Rwtot - total( random[ind].wsum[bin], /double)

          ;; rand
          sampval[j, bin] = Lsmod/Lwmod - Rsmod/Rwmod
      ENDFOR 

      jval[bin] = mean_check(sampval[*, bin], /double)

      lumdens[bin] = smean + (Nlenses-1)*(smean - jval[bin])
      
  ENDFOR 

  ;; Now jackknife covariance between variables
  FOR jbin=0L, Nbin-1 DO BEGIN 
      FOR ibin=jbin, Nbin-1 DO BEGIN 

          tmp = total( (sampval[*,ibin]-jval[ibin])*(sampval[*,jbin]-jval[jbin]), /double)
          covariance[jbin, ibin] = factor*tmp

          IF jbin NE ibin THEN covariance[ibin, jbin] = covariance[jbin, ibin]
          IF jbin EQ ibin THEN lumdenserr[ibin] = sqrt( covariance[ibin, ibin] )

      ENDFOR 
  ENDFOR 

  sampval=0
  jval=0

return
END 

