PRO e1e2boot, e1, e2, Nresamp, $
              me1, me2, me1err, me2err,$
              me1e1err, me1e2err, me2e2err, $
              me1e1err_err, me1e2err_err, me2e2err_err, $
              st_dev=st_dev, cut=cut

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: e1e2boot, e1, e2, Nresamp,me1, me2, me1err, me2err,me1e1err, me1e2err, me2e2err,me1e1err_err, me1e2err_err, me2e2err_err, st_dev=st_dev, cut=cut'
      print
      print,' If /sdev, also measure the me1e1err stuff'
      return
  ENDIF 

  ;; for simulations to test errors, so we want the variance
  ;; which gives the error on individual measurements

  ;; me1err is how well we determined mean
  ;; me1e1err is the sdev
  ;; me1e1err_err is the error in the sdev


  Nmeas = n_elements(e1)
  Nvar = 2

  IF keyword_set(cut) THEN BEGIN 
      w1=where( e1 EQ e1 AND e2 EQ e2)

      w=where( (e1[w1] LT 1.) AND (e1[w1] GT -1.) AND $
               (e2[w1] LT 1.) AND (e2[w1] GT -1.), Nmeas)
      IF Nmeas EQ 0 THEN BEGIN 
          print,'No e1/e2 passed cuts'
          return
      ENDIF 
      w=w1[w]
  ENDIF ELSE w=lindgen(Nmeas)

  ;; bootstrap me1, me2 sdev's
  data = fltarr(Nmeas, Nvar)
  data[*,0] = e1[w]
  data[*,1] = e2[w]

  print
  print,'Bootstrapping for mean e1,e2'
  print,'Creating sample'
  tt=systime(1)
  resample_data, data, Nmeas, Nvar, Nresamp, resamp
  ptime,systime(1)-tt
  print,'Bootstrapping'
  tt=systime(1)
  bootstrap, resamp, Nresamp, Nvar, meane, meaneerr
  ptime,systime(1)-tt

  me1 = meane[0]
  me2 = meane[1]
  me1err = meaneerr[0]
  me2err = meaneerr[1]

  IF NOT keyword_set(st_dev) THEN return
  print,'Bootstrapping for covariance(e1,e2)'
  print,'Creating sample'
  tt=systime(1)
  resample_data_cov, data, Nmeas, Nvar, Nresamp, resamp_cov, /st_dev
  ptime,systime(1)-tt
  print,'Bootstrapping'
  tt=systime(1)
  bootstrap, resamp_cov, Nresamp, Nvar*Nvar, tcov, tcoverr, tcovcov
  ptime,systime(1)-tt

  reconstruct_resamp_cov, tcov, tcoverr, Nvar, cov, coverr

  me1e1err = sqrt(cov[0,0])
  me1e2err = sqrt(abs(cov[0,1]))*sign(cov[0,1])
  me2e2err = sqrt(cov[0,0])

  me1e1err_err = coverr[0,0]/2./me1e1err
  me1e2err_err = coverr[0,1]/2./me1e2err
  me2e2err_err = coverr[1,1]/2./me2e2err

END 
