PRO cjackknife, data, Nmeas, Nvar, datamean, dataerr, covariance, statistic=statistic

  IF n_params() LT 4 THEN BEGIN
      print,'-Syntax: cjackknife, data, Nmeas, Nvar, datamean, dataerr, covariance, statistic=statistic'
      print
      print,'statistic=1 for mean'
      print,'statistic=2 for sdev'
      print,'statistic=3 for variance'
      return
  ENDIF 

  Nmeas = long(Nmeas)
  Nvar  = long(Nvar)
  IF n_elements(statistic) EQ 0 THEN statistic = 1L $
  ELSE statistic = long(statistic)

  ;; returned means and covariance
  datamean = fltarr(Nvar)
  covariance = fltarr(Nvar, Nvar)

  ;; we will get from diagonal elements
  dataerr = fltarr(Nvar)

  entry = esheldon_config("JACK_ENTRY")
  sofile = esheldon_config("JACK_SOFILE")

  retval = call_external(value=[0b, 0b, 0b, 0b, 0b, 0b], $
                         sofile, entry,$
                         data, Nmeas, Nvar, statistic, $
                         datamean, covariance)

  FOR i=0L, Nvar-1 DO dataerr[i] = sqrt(covariance[i,i])

  IF Nvar EQ 1 THEN BEGIN 
      datamean = datamean[0]
      dataerr = dataerr[0]
  ENDIF 

END 
