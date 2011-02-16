PRO wjackknife, dataSum, wSum, mean, error, covariance

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: wjackknife, dataSum, wSum, mean, error, covariance'
      return
  ENDIF 

  data_size = size(dataSum)

  nSub = data_size[1]  
  nBin = data_size[2]

  mean = replicate(dataSum[0,0], nBin)
  mean[*] = 0.0

  error = mean
  covariance = replicate(dataSum[0,0], nBin, nBin)
  covariance[*] = 0.0

  ;; temporary variables
  wSumTot = mean
  dataSumTot = mean
  dataMean = mean

  dataSumMod = dataSum
  wSumMod = wSum
  modMean = dataSum

  meanOverMods = mean

  nSubUsed = lonarr(nBin)

  ;; Get overall sums and the mod array
  FOR bin=0L, nBin-1 DO BEGIN 

      w = where(dataSum[*, bin] NE 0, nSubBin)
      nsubUsed[bin] = nSubBin
;      IF nSubBin NE nSub THEN $
;        message,'Some samples are empty: bin='+ntostr(bin)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Sums over all sub-samples
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      dataSumTot[bin] = total(dataSum[*, bin], /double)
      wSumTot[bin] = total(wSum[*, bin], /double)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; The mean over all samples
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF wSumTot[bin] GT 0 THEN dataMean[bin] = dataSumTot[bin]/wSumTot[bin] 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; The mod array for this bin
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
      dataSumMod[*, bin] = dataSumTot[bin] - dataSum[*, bin]
      wSumMod[*, bin]    = wSumTot[bin]    - wSum[*, bin]

      wz = where(wSumMod[*,bin] GT 0,nwz)
      IF nwz NE 0 THEN BEGIN 
          modMean[wz, bin] = dataSumMod[wz, bin]/wSumMod[wz, bin]
      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Now the mean over the mod samples
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      meanOverMods[bin] = mean_check(modMean[*, bin], /double)
      mean[bin] = dataMean[bin] + (nSub-1)*(dataMean[bin] - meanOverMods[bin])

  ENDFOR 

  factor = (nSub - 1.0)/nSub
  tfactor = (nSubUsed - 1.0)/nSubUsed
  print
  print,'Factor: ',factor
  print,'Actual factor: ',tfactor

  FOR jBin=0L, nBin-1 DO BEGIN 
      
      FOR iBin=jBin, nBin-1 DO BEGIN 
          
          tmp = total( (modMean[*, iBin] - meanOverMods[iBin])*(modMean[*, jBin] - meanOverMods[jBin]), /double)
          
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
