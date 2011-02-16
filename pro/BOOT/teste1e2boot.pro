PRO teste1e2boot, struct, statistic

  IF n_elements(struct) EQ 0 THEN BEGIN 
      struct = mrdfits('test_ellip_errors.fit',1)
  ENDIF 

  posangle = 60.
  s2n = 10.
  ellip = 0.3
  w=where(struct.posangle eq posangle and $
          struct.s2n eq s2n AND $
          struct.e EQ ellip,nw)

  me1 = struct[w].me1
  me2 = struct[w].me2

  w1=where(me1 EQ me1 AND me2 EQ me2)

  w2=where(me1[w1] GT -1. AND me1[w1] LT 1. AND $
           me2[w1] GT -1. AND me2[w1] LT 1.,Nmeas)
  w2=w1[w2]

  data = fltarr(Nmeas, 2)
  data[*,0] = me1[w2[0:Nmeas-1]]
  data[*,1] = me2[w2[0:Nmeas-1]]

  nVar = 2L

  ;; C jackknife code
  print
  print,"Jackknifing Ccode"
  tt=systime(1)
  cjackknife, data, Nmeas, Nvar, jdatamean, jdataerr, jcovariance, $
              statistic=statistic
  ptime,systime(1)-tt,6

  ;; C jackknife code for sdev
  print
  print,"Jackknifing Ccode for sdev"
  tt=systime(1)
  cjackknife, data, Nmeas, Nvar, jsdevmean, jsdeverr, jsdevcovariance, $
              statistic=2
  ptime,systime(1)-tt,6

  ;; C jackknife code for var
  ;;print
  ;;print,"Jackknifing Ccode for variance"
  ;;tt=systime(1)
  ;;cjackknife, data, Nmeas, Nvar, jvarmean, jvarerr, jvarcovariance, $
  ;;            statistic=3
  ;;ptime,systime(1)-tt,6

  ;; IDL jackknife code
  print
  print,"Jackknifing IDL code"
  tt=systime(1)
  jackknife, data, Nmeas, Nvar, jdatamean2, jdataerr2, jcovariance2
  ptime,systime(1)-tt,5

  ;; bootstrap
  Nresamp=500L

  ;; C bootstrap code
;  print
;  print,"Bootstrapping C code"
;  tt=systime(1)
;  cbootstrap, data, Nmeas, Nvar, Nresamp, $
;              bdatamean, bdataerr, bcovariance, bsampval, statistic=statistic
;  ptime,systime(1)-tt,5

  ;; IDL bootstrap code
;  print
;  print,'Bootstrapping IDL code'
;  tt=systime(1)
;  resample_data, data, Nmeas, Nvar, Nresamp, bsampval
;  bootstrap, bsampval, Nresamp, Nvar, bdatamean2, bdataerr2, bcovariance2
;  ptime,systime(1)-tt,5

  print
  print,'Means'
  print,' Input:   '
  print,[struct[w].e1,struct[w].e2]
;  print,' C Boot:  '
;  print,bdatamean,bdataerr
;  print,' IDL Boot:'
;  print,bdatamean2,bdataerr2
  print,' C Jack:  '
  print,jdatamean,jdataerr
  print,' IDL Jack: '
  print,jdatamean2,jdataerr2
  print
  print,[struct[w].e1e1err, struct[w].e2e2err]
  print,jsdevmean
  print,jsdeverr


END 
