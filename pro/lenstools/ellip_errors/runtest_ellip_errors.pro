PRO runtest_ellip_errors, niter, posangles, ellip, s2n, struct

;  posangles = [-90.0, -75.0, -60.0, -45.0, -30.0, -15.0, $
;               0.0, $
;                15.0,  30.0,  45.0,  60.0,  75.0,  90.0]

;  ellip = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
  
;  s2n = [20.0, 40.0, 80.0, 100.0, 120.0, 140.0, 160.0, 170.0, 180.0, 200.0]

  npos = n_elements(posangles)
  nellip = n_elements(ellip)
  ns2n = n_elements(s2n)

  arrval = fltarr(niter)
  ntot = npos*nellip*ns2n
  tstruct = create_struct('e', 0.0,$ ;input values
                          'e1', 0.0, $
                          'e2', 0.0, $
                          'Ixx',0.0,$
                          'Iyy',0.0,$
                          'Ixy',0.0,$
                          'e1e1err',0.0,$
                          'e1e2err',0.0,$
                          'e2e2err',0.0,$
                          'posangle', 0.0, $
                          's2n', 0.0, $
                          'T', 0.0, $
                          'niter', 0L, $
                          'nuse', 0L,$
                          $     ;measured values
                          'me1', arrval, $
                          'me2', arrval, $
                          'mT',  arrval,$
                          $     ;errors on measured values
                          'm_e1e1err', 0.0,$
                          'm_e1e2err', 0.0,$
                          'm_e2e2err', 0.0,$
                          'm_e1e1err_err',0.0,$
                          'm_e1e2err_err',0.0,$
                          'm_e2e2err_err',0.0,$
                          $     ;formal errors
                          'f_e1e1err', arrval,$
                          'f_e1e2err', arrval,$
                          'f_e2e2err', arrval,$
                          $     ;mean on formal errors
                          'f_e1e1err_mean',0.0,$
                          'f_e1e2err_mean',0.0,$
                          'f_e2e2err_mean',0.0,$
                          'f_e1e1err_err',0.0,$
                          'f_e1e2err_err',0.0,$
                          'f_e2e2err_err',0.0)

  struct = replicate(tstruct, ntot)

  ;; object sizes: ixx+iyy
  T = 10.0

  tt=systime(1)

  itot = 0L
  FOR ipos = 0L, npos-1 DO BEGIN 
      posangle = posangles[ipos]
      FOR iell = 0L, nellip-1 DO BEGIN 
          e = ellip[iell]
          FOR is2n = 0L, ns2n-1 DO BEGIN 
              
              ss2n = s2n[is2n]

              tt=systime(1)
              test_ellip_errors, e, posangle, ss2n, T, niter, $
                                 e1, e2, Ixx, Iyy, Ixy, $
                                 e1e1err, e1e2err,e2e2err,Terr,$
                                 me1, me2, mT, $
                                 f_e1e1err, f_e1e2err, f_e2e2err, fTerr


              ;; failures -9999.
              w1=where(me1 EQ me1 AND me2 EQ me2)

              w=where(me1[w1] GT -1. AND me1[w1] LT 1. AND $
                      me2[w1] GT -1. AND me2[w1] LT 1.,nuse)
              w=w1[w]
              fac = 1./(nuse-1.)

              e1var = fac*total( (me1[w]-e1)^2 )
              e2var = fac*total( (me2[w]-e2)^2 )
              e1e2var = fac*total( (me1[w]-e1)*(me2[w]-e2) )
              Tvar = fac*total( (mT[w]-T)^2 )

              ;; use bootstrap to get the covariances
              ;; and errors on those covariances

              nresamp=1000L
              e1e2boot, me1[w], me2[w], nresamp, $
                        tmpe1, tmpe2, tmpe1err, tmpe2err,$
                        m_e1e1err, m_e1e2err, m_e2e2err, $
                        m_e1e1err_err, m_e1e2err_err, m_e2e2err_err,/st_dev

              ;; average and err on avg. for the
              ;; formal errors
              f_e1e1err_mean = mean( f_e1e1err[w] )
              f_e1e2err_mean = mean( f_e1e2err[w] )
              f_e2e2err_mean = mean( f_e2e2err[w] )

              f_e1e1err_err = sdev( f_e1e1err[w] )
              f_e1e2err_err = sdev( f_e1e2err[w] )
              f_e2e2err_err = sdev( f_e2e2err[w] )

              ;; The input values
              struct[itot].e = e
              struct[itot].e1 = e1
              struct[itot].e2 = e2
              struct[itot].posangle = posangle
              struct[itot].s2n = ss2n
              struct[itot].T = T
              struct[itot].niter = niter
              struct[itot].nuse = nuse
              struct[itot].e1e1err = e1e1err
              struct[itot].e1e2err = e1e2err
              struct[itot].e2e2err = e2e2err

              ;; measured e's and the error
              struct[itot].me1 = me1
              struct[itot].me2 = me2
              struct[itot].mT = mT
              struct[itot].m_e1e1err = m_e1e1err
              struct[itot].m_e1e2err = m_e1e2err
              struct[itot].m_e2e2err = m_e2e2err
              struct[itot].m_e1e1err_err = m_e1e1err_err
              struct[itot].m_e1e2err_err = m_e1e2err_err
              struct[itot].m_e2e2err_err = m_e2e2err_err

              ;; formal errors and means
              struct[itot].f_e1e1err = f_e1e1err
              struct[itot].f_e1e2err = f_e1e2err
              struct[itot].f_e2e2err = f_e2e2err
              struct[itot].f_e1e1err_mean = f_e1e1err_mean
              struct[itot].f_e1e2err_mean = f_e1e2err_mean
              struct[itot].f_e2e2err_mean = f_e2e2err_mean
              struct[itot].f_e1e1err_err = f_e1e1err_err
              struct[itot].f_e1e2err_err = f_e1e2err_err
              struct[itot].f_e2e2err_err = f_e2e2err_err


              itot = itot+1L

              ptime,systime(1)-tt

          ENDFOR 
      ENDFOR 
  ENDFOR 

  ptime, systime(1)-tt


END 
