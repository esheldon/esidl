PRO resample_data_cov_getmean, array, mean, weights=weights

  IF n_elements(weights) NE 0 THEN BEGIN 

      wsum = total(weights)
      mean = total(array*weights)/wsum

  ENDIF ELSE BEGIN 

      mean = mean_check(array)

  ENDELSE 

END 

PRO resample_data_cov_getcov, array1, mean1, array2, mean2, cov, $
                              weights1=weights1, weights2=weights2

  narray1 = n_elements(array1)
  narray2 = n_elements(array2)
  IF narray1 NE narray2 THEN message,'array1 and array2 must be same size'
  IF n_elements(weights1) NE 0 THEN BEGIN 

      weights = sqrt(weights1*weights2)
      wsum = total(weights)
      cov = total(weights*(array1-mean1)*(array2-mean2))/wsum

  ENDIF ELSE BEGIN 

      cov = total( (array1-mean1)*(array2-mean2) )/(narray1-1)

  ENDELSE 

END 

PRO resample_data_cov, data, Nmeas, Nvar, nresamp, $
                       resamp_cov, weights=weights, $
                       st_dev=st_dev

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: resample_data, data, Nmeas, Nvar, nresamp, resamp_cov, weights=weights'
      print,' data must be of the form array(Nmeasurements,Nvariables)'
      print,' weights should be same size as data'
      print,' the mean over the resamples'
      print
      print,' The output can be sent to bootstrap.pro'
      print,' Use reconstruct_resamp_cov.pro to get covariance matrix'
      return
  ENDIF 

  ;; data must be of the form array[N,M]
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

  doweight=0
  IF n_elements(weights) NE 0 THEN doweight = 1
      
  ;;Nunique = (Nvar^2 + Nvar)/2 = Nvar*(Nvar + 1)/2
  resamp_cov = replicate(data[0], nresamp, Nvar*Nvar)

  FOR samp = 0L, nresamp-1 DO BEGIN 

      ;; get a random set of objects with replacement
      sind = round( (Nmeas-1)*randomu(seed, Nmeas) )

      FOR jvar=0L, Nvar-1 DO BEGIN 

          IF doweight THEN BEGIN
              resample_data_cov_getmean,data[sind,jvar],jmean,$
                                        weights=weights[sind,jvar]
          ENDIF ELSE BEGIN 
              resample_data_cov_getmean,data[sind,jvar],jmean
          ENDELSE 

          FOR ivar=jvar, Nvar-1 DO BEGIN 

              IF ivar EQ jvar THEN imean = jmean ELSE BEGIN 
                  IF doweight THEN BEGIN
                      resample_data_cov_getmean,data[sind,ivar],imean,$
                                                weights=weights[sind,ivar]
                  ENDIF ELSE BEGIN 
                      resample_data_cov_getmean,data[sind,ivar],imean
                  ENDELSE 
              ENDELSE 
              IF doweight THEN BEGIN
 
                  resample_data_cov_getcov,data[sind,jvar],jmean,$
                                           data[sind,ivar],imean,$
                                           tcov, $
                                           weights1=weights[sind,jvar],$
                                           weights2=weights[sind,ivar]
              ENDIF ELSE BEGIN 
                  resample_data_cov_getcov,data[sind,jvar],jmean,$
                                           data[sind,ivar],imean,$
                                           tcov
              ENDELSE 
              ;; like filling resamp[ivar,jvar] if was 2-D form
              IF keyword_set(st_dev) THEN BEGIN
                  cov=sqrt(abs(tcov))*sign(tcov)
              ENDIF ELSE BEGIN 
                  cov=tcov
              ENDELSE 
              resamp_cov[samp, jvar*Nvar + ivar] = cov
              IF ivar NE jvar THEN resamp_cov[samp, ivar*Nvar + jvar] = cov

          ENDFOR 
      ENDFOR 
  ENDFOR 


END 
