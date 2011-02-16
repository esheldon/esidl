PRO test_e1e2err, e, posangle, s2n, T, e1, e2, Ixx,Iyy,Ixy,e1e1err,e1e2err,e2e2err,Terr,te1e1err,te1e2err,te2e2err,tTerr

  ;; see what happens to e1e2err when add noise

  ntest = 10000L

  defval=-9999d
  te1e1err=replicate(defval,ntest)
  te1e2err=replicate(defval,ntest)
  te2e2err=replicate(defval,ntest)
  tTerr=replicate(defval,ntest)

  theta = posangle*!dpi/180d

  e1 = e*cos(2.*theta)
  e2 = e*sin(2.*theta)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; make object along x-axis with given ellipticity
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  
  ;; ensure iyy = aratio^2*ixx
  Ixx = (T/2.)*(1 + e1)
  Iyy = (T/2.)*(1 - e1)
  Ixy = (T/2.)*e2

  det = Ixx*Iyy - Ixy^2

  ;; use sky = 100.
  sky = 100.

  ;; make the image without noise
  nx = 25L
  ny = 25L

  l  = lindgen(nx*ny)
  x  = reform( float(l mod nx), nx,ny)
  y  = reform( float(l / nx),   nx,ny) 
  cenx = (nx-1.)/2.
  ceny = (ny-1.)/2.

  ;; from s2n: Dave's code uses all the pixels!
  ;; only use skyh noise
;  counts = s2n*sqrt( sky*nx*ny )
  counts = s2n*sqrt( sky*!pi*16.0*sqrt(det) )
  Amp = counts/2./!pi/sqrt(det)

  ;; get predicted errors
  daves_mom_errors,Amp,Ixx,Iyy,Ixy,sky,covmat,tVar_e1,tVar_e2,tVar_T,tCov_e1e2,tCorr
  e1e1err=sqrt(tVar_e1)
  e1e2err=sqrt(abs(tCov_e1e2))*sign(tCov_e1e2)
  e2e2err=sqrt(tVar_e2)
  Terr=sqrt(tVar_T)

  ;; now add noise at 10% level to ixx, iyy, ixy
  errsize = 0.2
  Ixx_err = (T/2.)*randomn(seed, ntest)*errsize
  Iyy_err = (T/2.)*randomn(seed, ntest)*errsize
  Ixy_err = (T/2.)*randomn(seed, ntest)*errsize

  FOR i=0L, ntest-1 DO BEGIN 
      
      tIxx = Ixx + Ixx_err[i]
      tIyy = Iyy + Iyy_err[i]
      tIxy = Ixy + Ixy_err[i]

      tDet = tIxx*tIyy - tIxy^2

      tAmp = counts/2./!pi/sqrt(tDet)

      daves_mom_errors,tAmp,tIxx,tIyy,tIxy,sky,covmat,$
                       Var_e1,Var_e2,Var_T,Cov_e1e2,Corr

      IF n_elements(Var_e1) NE 0 THEN BEGIN 
          te1e1err[i] = sqrt(Var_e1)
          te1e2err[i] = sqrt(abs(Cov_e1e2))*sign(Cov_e1e2)
          te2e2err[i] = sqrt(Var_e2)
          tTerr[i]    = sqrt(Var_T)
      ENDIF 

  ENDFOR 

  w=where(te1e2err GT defval)

  pcharold=!p.charsize
  !p.multi=[0,0,3]
  !p.charsize=2.0

  bin=0.0001
  plothist,te1e1err[w],bin=bin,xtitle='e1e1err'
  oplot,[e1e1err,e1e1err],[0,10000],color=!red,thick=5
  plothist,te1e2err[w],bin=bin,xtitle='e1e2err'
  oplot,[e1e2err,e1e2err],[0,10000],color=!red,thick=5
  plothist,te2e2err[w],bin=bin,xtitle='e2e2err'
  oplot,[e2e2err,e2e2err],[0,10000],color=!red,thick=5

  !p.charsize=pcharold

END 
