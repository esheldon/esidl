PRO wbootstrap, dataSum, wSum, nResample, mean, error, covariance, meanSample=meanSample

  IF n_params() LT 6 THEN BEGIN 
      print,'-Syntax: wbootstrap, dataSum, wSum, nResample, mean, error, covariance'
      return
  ENDIF 

  data_size = size(dataSum)

  nSub = data_size[1]  
  nBin = data_size[2]

  mean = dblarr(nBin)
  mean[*] = 0.0

  error = mean
  covariance = dblarr(nBin, nBin)
  covariance[*] = 0.0

  ;; temporary variables
  wSumTot = mean
  dataSumTot = mean
  dataMean = mean

  dataSumResample = dblarr(nResample, nBin)
  wSumResample = dataSumResample
  meanSample = dataSumResample
  meanOverSamples = dblarr(nBin)

  nSubUsed = lonarr(nBin)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; create samples and the means in
  ;; each bin for those samples
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR iSamp=0L, nResample-1 DO BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Draw with replacement nResample times
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      sind = round( (nSub-1)*randomu(seed, nSub) )

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Means in each bin
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;

      FOR bin=0L, nBin-1 DO BEGIN 

          dataSumResample[iSamp, bin] = total(dataSum[sind, bin], /double)
          wSumResample[iSamp, bin]    = total(wSum[sind, bin], /double)

          meanSample[iSamp, bin] = dataSumResample[isamp, bin]/wSumResample[iSamp, bin]

      ENDFOR 

  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Get overall sums and means
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  

  FOR bin=0L, nBin-1 DO BEGIN 

      w = where(dataSum[*, bin] NE 0, nSubBin)
      nSubUsed[bin] = nSubBin

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Sums over all measurements
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      dataSumTot[bin] = total(dataSum[*, bin], /double)
      wSumTot[bin]    = total(wSum[*, bin], /double)
      dataMean[bin] = dataSumTot[bin]/wSumTot[bin] 

      meanOverSamples[bin] = $
        total( dataSumResample[*, bin], /double ) / total( wSumResample[*, bin], /double )

  ENDFOR 
  
  factor = 1./(nResample-1.)
  FOR jBin=0L, nBin-1 DO BEGIN 
      
      FOR iBin=jBin, nBin-1 DO BEGIN 
          
          tmp = total( (meanSample[*, iBin] - meanOverSamples[iBin])*(meanSample[*, jBin] - meanOverSamples[jBin]) , /double)
          
          covariance[jBin, iBin] = factor*tmp
          
          IF jBin NE iBin THEN BEGIN 
              covariance[iBin, jBin] = covariance[jBin, iBin]
          ENDIF 
          IF jBin EQ iBin THEN BEGIN 
              error[iBin] = sqrt(covariance[iBin, iBin])
          ENDIF 
          
      ENDFOR 
      
  ENDFOR 

END 
