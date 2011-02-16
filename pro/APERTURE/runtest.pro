PRO runtest, n

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: runtest, n'
      print,' n is number for procedure m_ap_test[n].pro'
      return
  ENDIF 

  ;; Read bhuv's unsmoothed and smoothed maps

  dir = '/sdss4/data1/esheldon/LARGEMAPS/'
  bhuvk=mrdfits(dir + 'gamma1_zs0304_cic.fits')
  bhuvksm = mrdfits(dir + 'gamma1_zs0304_cicsm.fits')

  minx = 0.                     ;degrees
  maxx = 5.
  miny = 0.
  maxy = 5.

  CASE n OF 
      1: BEGIN
          m_ap_test1, bhuvk, minx, maxx, miny, maxy, $
                      m_ap2, rsize, M_Power, lind, m_ap2err, M_Powererr
          m_ap_test1, bhuvksm, minx, maxx, miny, maxy, $
                      sm_ap2, srsize, sM_Power, slind, sm_ap2err, sM_Powererr
      END 
      2: BEGIN
          m_ap_test2, bhuvk, minx, maxx, miny, maxy, $
                      m_ap2, rsize, M_Power, lind, m_ap2err, M_Powererr
          m_ap_test2, bhuvksm, minx, maxx, miny, maxy, $
                      sm_ap2, srsize, sM_Power, slind, sm_ap2err, sM_Powererr
      END 
      3: BEGIN 
          m_ap_test3, bhuvk, minx, maxx, miny, maxy, $
                      m_ap2, rsize, M_Power, lind, m_ap2err, M_Powererr
          m_ap_test3, bhuvksm, minx, maxx, miny, maxy, $
                      sm_ap2, srsize, sM_Power, slind, sm_ap2err, sM_Powererr
      END 
      ELSE : BEGIN
          print,ntostr(n),' is not a valid  number'
          return
      END 
  ENDCASE 

  title = 'Simulated Kappa Maps'
  xtitle1 = 'Arcminutes'
  xtitle2 = 'l    (= 2*pi/theta  theta in radians)'
  ytitle1 = '<M_ap^2>'
  ytitle2 = 'P(kappa)'

  message = ['Raw', 'Smoothed']
  psym = [1, 0]

  ploterr, rsize, m_ap2, m_ap2err, /xlog, psym=1, $
        title=title, xtitle=xtitle1, ytitle=ytitle1
  oploterr, srsize, sm_ap2, sm_ap2err
  legend, message, psym=psym

  key=get_kbrd(1)

  ploterr, lind, M_Power, M_Powererr, /xlog, /ylog, psym=2, $
        title=title, xtitle=xtitle2, ytitle=ytitle2
  oploterr, slind, sM_Power, sM_Powererr
  legend, message, psym=psym

  return 
END 
