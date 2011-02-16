PRO cbootstrap, data, Nmeas, Nvar, Nresamp, datamean, dataerr, covariance, sampval, statistic=statistic

  IF n_params() LT 5 THEN BEGIN
      print,'-Syntax: cbootstrap, data, Nmeas, Nvar, statistic, Nresamp, datamean, dataerr, covariance, sampval, statistic=statistic'
      return
  ENDIF 

  Nmeas = long(Nmeas)
  Nvar  = long(Nvar)
  IF n_elements(statistic) EQ 0 THEN statistic = 1L $
  ELSE statistic = long(statistic)
  Nresamp = long(Nresamp)

  ;; returned means and covariance
  datamean = fltarr(Nvar)
  covariance = fltarr(Nvar, Nvar)

  ;; we will get from diagonal elements
  dataerr = fltarr(Nvar)

  ;; returned statistic from random samples
  sampval = fltarr(Nresamp, Nvar)

  sofile = '/net/cheops2/home/esheldon/ccode/bootstrap/call_bootstrap.so'
  entry = 'call_bootstrap'

  retval = call_external(value=[0b, 0b, 0b, 0b, 0b, 0b, 0b, 0b], $
                         sofile, entry,$
                         data, Nmeas, Nvar, statistic, $
                         sampval, Nresamp, $
                         datamean, covariance)

  FOR i=0L, Nvar-1 DO dataerr[i] = sqrt(covariance[i,i])

END 
