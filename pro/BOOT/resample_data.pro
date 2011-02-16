PRO resample_data, data, Nmeas, Nvar, nresamp, $
                   resamp, weights=weights, sdev=sdev, sums=sums

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: resample_data, data, Nmeas, Nvar, nresamp, resamp, weights=weights, sdev=sdev'
      print,'data must be of the form array(Nmeasurements,Nvariables)'
      print,'weights should be same size as data'
      print,' If /sdev, then resamp contains the standard dev. rather than'
      print,' the mean over the resamples'
      print
      print,'The output can be sent to bootstrap.pro'
      return
  ENDIF 

  ;; data must be of the form array[M,N]
  ;; where N is the number of variables and M is the
  ;; number of measurements.
  ;; 
  ;; e.g. Lets say we have Ngal galaxies and
  ;; we measure the number of neighbors in Nbin radial bins.  THen
  ;; the data array is an array of "number of neighbors" which is
  ;; of size (Nbin, Ngal)
  ;;
  ;; Then we resample:  nresamp times we choose Ngal objects randomly
  ;; with replacement.  The mean in each bin is then calculated, possibly
  ;; with input weights.  The output is then an array (Nresamp, Nbin)
  ;;
  ;; another example: say we did a simulation and measured e1,e2 N times
  ;; for N different noise realizations.  Then Nmeas = N and Nvar =2

  ;; Note: can send /sdev and it will do standard deviation rather
  ;; than mean, but to get covariances you should use resample_data_cov

  IF n_elements(weights) NE 0 THEN doweight=1 ELSE doweight=0
  
  resamp = replicate(data[0], nresamp, nvar)

  FOR samp = 0L, nresamp-1 DO BEGIN 

      ;; get a random set of Measurements with replacement
      sind = round( (Nmeas-1)*randomu(seed, Nmeas) )

      FOR var=0L, nvar-1 DO BEGIN 

          IF doweight THEN BEGIN
 
              ;; weighted average
              wsum = total(weights[sind,var])
              mean = total(data[sind,var]*weights[sind,var])/wsum

              IF NOT keyword_set(sdev) THEN BEGIN 
                  resamp[samp, var] = mean
              ENDIF ELSE BEGIN
                  tVar = total( weights[sind,var]*(data[sind,var]-mean)^2 )
                  tVar = tVar/wsum
                  resamp[samp,var] = sqrt(tVar)

              ENDELSE 
          ENDIF ELSE BEGIN 

              ;; unweighted average
              mean = mean_check( data[sind,var] )
              
              IF NOT keyword_set(sdev) THEN BEGIN 
                  resamp[samp, var] = mean
              ENDIF ELSE BEGIN 
                  tVar = total( (data[sind,var] - mean)^2 )/(Nmeas-1.)
                  resamp[samp, var] = sqrt(tVar)
              ENDELSE 

          ENDELSE 
      ENDFOR 
  ENDFOR 


END 
