PRO jackknife, data, datamean, dataerr, covariance

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: jackknife, data, datamean, dataerr, covariance'
      print
      return
  ENDIF 

  tt=systime(1)

  ;; See page Section 6.5 Lupton
  ;;
  ;; Note this does mean but you could make it
  ;; some other statistic as well and get
  ;; variance on it: e.g. variance and variance error
  ;;
  ;; No Weights Yet!!

  sz = size(data)
  IF sz[0] NE 2 THEN message,'data must be a NmeasXNvar array. E.g. NmeasXNradius'

  Nmeas = sz[1]
  Nvar  = sz[2]

  ;; Outputs: mean, covariance and diagonal terms
  datamean = dblarr(Nvar)
  dataerr = datamean
  covariance = dblarr(Nvar, Nvar)

  factor = (Nmeas - 1.)/Nmeas

  sampval = dblarr(Nmeas, Nvar)
  jval = dblarr(Nvar)
  stot = 0d

  ;; measure totals. Just subtract off appropriate thing
  ;; to get means in subsamples

  FOR ivar=0L, Nvar-1 DO BEGIN

      stot = total(data[*, ivar], /double)
      smean = stot/Nmeas

      FOR j=0L, Nmeas-1 DO BEGIN 
          sampval[j, ivar] = ( stot - data[j, ivar] )/(Nmeas-1.)
      ENDFOR 

      jval[ivar] = mean_check(sampval[*, ivar],/double)

      ;; this is bias-corrected valuex
      datamean[ivar] = smean + (Nmeas-1)*(smean - jval[ivar])
  ENDFOR 

  ;; Now jackknife covariance between variables
  FOR jvar=0L, Nvar-1 DO BEGIN 
      FOR ivar=jvar, Nvar-1 DO BEGIN 

          tmp = total( (sampval[*,ivar]-jval[ivar])*(sampval[*,jvar]-jval[jvar]), /double)
          covariance[jvar, ivar] = factor*tmp

          IF jvar NE ivar THEN covariance[ivar, jvar] = covariance[jvar, ivar]
          IF jvar EQ ivar THEN dataerr[ivar] = sqrt( covariance[ivar, ivar] )

      ENDFOR 
  ENDFOR 
  
  sampval=0
  jval=0

return

END 
