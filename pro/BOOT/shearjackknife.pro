PRO shearjackknife, lenses, random, sigma, sigmaerr, covariance, wuse=wuse

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: shearjackknife, lenses, random, sigma, sigmaerr, covariance, wuse=wuse'
      return
  ENDIF 

  Nlenses = n_elements(lenses)
  Nrand = n_elements(random)
  IF n_elements(wuse) NE 0 THEN BEGIN 
      Nbin = n_elements(wuse)
  ENDIF ELSE BEGIN 
      Nbin = n_elements(lenses[0].rsum)
      wuse = lindgen(nbin)
  ENDELSE  

  sigma  = dblarr(Nbin)
  sigmaerr = sigma
  covariance = dblarr(Nbin, Nbin)

  factor = (Nlenses - 1.)/Nlenses

  sampval = dblarr(Nlenses, Nbin)
  jval = dblarr(Nbin)
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
  print,'Number of random: ',Nrand
  print,'Jackknifing'

  ;; Use histogram to get random for each lens
  ;; using histogram this way, with min=0, forces a bin for every
  ;; integer from 0 to the max. That way, we can use
  ;; arr1[i] as subscript for reverse_indices! 
  
  ;; match multi gets all matches, but also returnes reverse_indices
  ;; for our use.  Assumes each lens zindex is represented in random
;  match_multi, lenses.zindex, random.zindex, rmatch, reverse_indices=revind
  
  ;; this gets all matches

  ssh = total(lenses.sshsum, /double)/total(lenses.wsum_ssh, /double)

  FOR bin=0L, Nbin-1 DO BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; correction factor for this radial bin
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      corr = total(lenses.wsum[bin], /double)/total(random.wsum[bin], /doublex)*float(Nrand)/nLenses/ssh

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; lenses which contributed to this bin
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      terr = reform(lenses.sigmaerr[bin])
      w=where(terr NE 0.0 AND finite(terr), Ngood)

      ;; matching random points
      match_multi, lenses[w].zindex, random.zindex, rmatch, reverse_indices=revind

      terr = terr[w]
      weights = 1./terr^2

      sigsum_tot = total( lenses[w].sigma[bin]*weights, /double)
      wtot = total( weights, /double)
      sigmean = sigsum_tot/wtot*corr

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; loop over the subsamples
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      FOR j=0L, Ngood-1 DO BEGIN 

          sigmod = sigsum_tot - lenses[w[j]].sigma[bin]*weights[j]
          wmod = wtot - weights[j]

          sampval[w[j], bin] = sigmod/wmod*corr

      ENDFOR 

      jval[bin] = mean_check(sampval[w, bin], /double)

      sigma[bin] = sigmean + (Ngood-1)*(sigmean - jval[bin])
      
  ENDFOR 

  ;; Now jackknife covariance between variables
  FOR jbin=0L, Nbin-1 DO BEGIN 
      wj=where(sampval[*,jbin] NE 0)
      FOR ibin=jbin, Nbin-1 DO BEGIN 

          wi=where(sampval[*,ibin] NE 0)

          tmp = total( (sampval[wi,ibin]-jval[ibin])*(sampval[wj,jbin]-jval[jbin]), /double)
          covariance[jbin, ibin] = factor*tmp

          IF jbin NE ibin THEN covariance[ibin, jbin] = covariance[jbin, ibin]
          IF jbin EQ ibin THEN sigmaerr[ibin] = sqrt( covariance[ibin, ibin] )

      ENDFOR 
  ENDFOR 

  sampval=0
  jval=0

return
END 

