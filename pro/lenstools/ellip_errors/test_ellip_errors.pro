FUNCTION test_ellip_gauss, x, y, p

  ;; evaluate a 2D gaussian
  smax=10D
  p=double(p)
  det=p(3)*p(4)-p(5)^2
  chisq= [p(4)*(x-p(0))^2-2.0*p(5)*(x-p(0))*(y-p(1))+p(3)*(y-p(1))^2]/(det> 1e-9) > 0.0
  mask = chisq LT (smax^2)
  ph=mask*exp(-.5* chisq *mask )

  ph = p[2]*ph/total(ph)

  return,ph

END 

PRO test_ellip_errors_amom, im0, sky, niter, e1, e2, T, e1e1err, e1e2err, e2e2err, Terr, imnoise=imnoise
  
  defval = -9999.
  e1 = replicate(defval,niter)
  e2 = e1
  T  = e1

  e1e1err=e1
  e1e2err=e1
  e2e2err=e1
  Terr=e1

  FOR i=0L, niter-1 DO BEGIN 

      IF keyword_set(imnoise) THEN BEGIN 
          im = im0+sky
          add_noise, im
      ENDIF ELSE BEGIN 
          sky_im = im0*0.0 + sky
          add_noise, sky_im
          im = im0 + sky_im
      ENDELSE 

      daves_ad_mom, im, sky, ixx, iyy, ixy, $
                    te1, te2, var_e1, var_e2, cov_e1e2, var_T
      
      IF n_elements(var_e1) NE 0 THEN BEGIN 
          e1e1err[i] = sqrt(var_e1)
          e1e2err[i] = sqrt(abs(cov_e1e2))*sign(cov_e1e2)
          e2e2err[i] = sqrt(var_e2)
          Terr[i] = sqrt(var_T)

          e1[i] = te1
          e2[i] = te2
          T[i] = ixx+iyy
      ENDIF 
      IF ( (i+1) MOD 1000) EQ 0 THEN print,ntostr(i+1)+'/'+ntostr(niter)

  ENDFOR 

  return

END 

PRO test_ellip_errors, e, posangle, s2n, T, niter, $
                       e1, e2, Ixx,Iyy,Ixy,e1e1err,e1e2err,e2e2err,Terr,$
                       me1, me2, mt, fe1e1err, fe1e2err, fe2e2err, fTerr,$
                       imnoise=imnoise

  ;; Test Dave's error formulas versus simulations

  ;; inputs:
  ;;  e: ellipticity
  ;;  posangle: in degrees
  ;;  s2n: signal to noise ratio
  ;;  niter: number of times to
  ;;         measure with noise

  ;; outputs:
  ;; e1, e2: the e1,e2 for input e and posangle
  ;; T: the T chosen below
  ;; me1: fltarr(niter) of measured values
  ;; me2:
  ;; mT:
  ;; e1e1err, etc., formal errors


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
  daves_mom_errors,Amp,Ixx,Iyy,Ixy,sky,covmat,Var_e1,Var_e2,Var_T,Cov_e1e2,Corr
  e1e1err=sqrt(Var_e1)
  e1e2err=sqrt(abs(Cov_e1e2))*sign(Cov_e1e2)
  e2e2err=sqrt(Var_e2)
  Terr=sqrt(Var_T)

  pars = [cenx, ceny, counts, Ixx, Iyy, Ixy]
  im0 = test_ellip_gauss(x, y, pars)

  print
  print,'Input e: '+ntostr(e)
  print,'Input e1: '+ntostr(e1)+' e2: '+ntostr(e2)
  print,'Position angle: '+ntostr(posangle)
  print,'S/N: '+ntostr(s2n)
  print,'sky: '+ntostr(sky)
  print,'counts: '+ntostr(counts)
  print,'Amplitude: '+ntostr(Amp)
  print,'e1e1err: '+ntostr(e1e1err)
  print,'e1e2err: '+ntostr(e1e2err)
  print,'e2e2err: '+ntostr(e2e2err)
  print,'Terr: '+ntostr(Terr)
  print

  !p.multi=[0,0,2]

  sky_im = im0*0.0 + sky
  add_noise, sky_im
  copyim = im0+sky_im
  tvim2, copyim, title='e: '+ntostr(e)+' PosAng: '+ntostr(posangle)+' S/N: '+ntostr(s2n)
  plot,im0
  !p.multi=0


  test_ellip_errors_amom, im0, sky, niter, me1, me2, mT, $
                          fe1e1err, fe1e2err, fe2e2err, fTerr, imnoise=imnoise

  

  return

END 
