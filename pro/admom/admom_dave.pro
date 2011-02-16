pro daves_ad_mom,im,sky,$
                 ixx,iyy,ixy,e1,e2,$
                 Var_e1,Var_e2,Cov_e1e2,Var_T,$
                 break=BREAK,amp=amp,noerr=noerr
  
  if n_params() eq 0 then begin
      print,'-syntax am_mom,im,p,s,sky,ixx,iyy,ixy,mx,my,c,im2,cov,err'
      return
  endif
  
  defval = -9999.

  maxiter=20
  siz=size(im)
  nx=long(siz(1))
  ny=long(siz(2))
  npix=nx*ny
  
  l=lindgen(npix)
  d=double(im(*))
  
  x=double(l mod nx)
  y=double(l /nx)
  
  c=total(d)
  s=(sky*npix)/c
  p=(1d)-s
  d=d-sky
  
  ixx=20d
  iyy=20d
  ixy=0d
  mx=(nx-1.)/(2d)
  my=(ny-1.)/(2d)
  
  det=ixx*iyy-ixy^2
  chi2= [iyy*(x-mx)^2-2.0*ixy*(x-mx)*(y-my)+ixx*(y-my)^2]/det
  
  phi=exp(-.5*(chi2 < 200d)  )
  phi=phi/total(phi)

  for i=0L, maxiter do begin
      det=ixx*iyy-ixy^2
      chi2= [iyy*(x-mx)^2-2.0*ixy*(x-mx)*(y-my)+ixx*(y-my)^2]/det
      phi=exp(-.5*(chi2 < 200d)  )
      phi=phi/total(phi)
      ;; w=p*phi/(p*phi+(s/npix))
      ;; r=sqrt((x-mx)^2+(y-my)^2)
      ;; plot,r,w,psym=3
      dw=d*phi
      p=total(dw)/total(phi^2)
      p=p/c
      ;; s=(1d)-p
      den=total(dw)
      mx=total(x*dw)/den
      my=total(y*dw)/den
      ixx=2.0*total(((x-mx)^2)*dw)/den
      ixy=2.0*total((x-mx)*(y-my)*dw)/den
      iyy=2.0*total(((y-my)^2)*dw)/den

  ENDFOR
  
  det=ixx*iyy-ixy^2
  m=p*c

  e1=(ixx-iyy)/(ixx+iyy)
  e2=2*ixy/(ixx+iyy)

  IF keyword_set(noerr) THEN return

  Amp=m/(2.0*!pi*sqrt(det))
  
  IF Amp GT 0.0 THEN BEGIN 
      ;;print,'amp',amp
      daves_mom_errors,Amp,Ixx,Iyy,Ixy,sky,covmat,Var_e1,Var_e2,Var_T,Cov_e1e2,Corr
      
  ENDIF ELSE BEGIN 

     return

  ENDELSE 
  return
end










