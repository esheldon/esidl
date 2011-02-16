PRO daves_mom_errors,A,Ixx,Iyy,Ixy,sky,cov,Var_e1,Var_e2,Var_T,Cov_e1_e2,Corr,$
print=print,det=det
;calculate the errors on the adaptive moments
;the sky variance is sky and
;Ixx,Iyy,Ixy are the obvious moments
;A is the amplitude (not total counts)
;this can just be cut and pasted into a C program

  if n_params() eq 0 then begin
      print,'-syntax mom_errors,A,Ixx,Iyy,Ixy,sky,cov,Var_e1,Var_e2,Var_T,Cov_e1_e2,Corr,print=print'
      return
  endif
  
  delvarx, cov, Var_e1, Var_e2, Var_T, cov_e1_e2, Corr

  pi=3.14159
  D=Ixx*Iyy-Ixy*Ixy
  if D le 0 then begin
      print,'Determinant negative, returning'
      return
      ;; should do something better
      ;; than just return but fine for now
  endif
                                ;a normalization factor
  F=pi*sqrt(D)/sky

                                ;now calculate the 10 independent elements 
                                ;of the 4x4 Fischer matrix 
  Fm=dblarr(4,4)	;IDL command that makes a 4x4 matrix (double precision)
  
  Fm(0,0)=F
  fac=F*A/(4.0*D)
  Fm(0,1)=fac*Iyy
  Fm(1,0)=Fm(0,1)
  Fm(0,2)=fac*Ixx	
  Fm(2,0)=Fm(0,2)
  Fm(0,3)=-fac*2.0*Ixy	
  Fm(3,0)=Fm(0,3)
  
  fac=3.0*F*A*A/(16.0*D*D)
  Fm(1,1)=fac*Iyy*Iyy
  Fm(2,2)=fac*Ixx*Ixx
  Fm(3,3)=fac*4.0*(Ixy*Ixy+D/3.0)
  
  Fm(1,2)=Fm(3,3)/4.0
  Fm(2,1)=Fm(1,2)
  Fm(1,3)=fac*(-2.0*Iyy*Ixy)
  Fm(3,1)=Fm(1,3)
  Fm(2,3)=fac*(-2.0*Ixx*Ixy)
  Fm(3,2)=Fm(2,3)
                                ;We have a matrix. Now invert it.
  
                                ;check to see if it is positive definite
                                ;it had better be
  det=determ(Fm)
  if det le 0 then begin
      print,'Fischer matrix not positive definite'
      ;;wait,5
      return
  endif
  
  cov=invert(Fm,/double)
;cov=Fm
                                ;IDL command that inverts a matrix
  
;print,'total',total(Fm)
  T=Ixx+Iyy
  e1=(Ixx-Iyy)/T
  e2=2.0*Ixy/T
  
  T4=4.0/(T*T*T*T)
;print,"T4 ",T4
  
;print,'here',(cov(1,1)*Iyy*Iyy +cov(2,2)*Ixx*Ixx -2.0*Ixx*Iyy*cov(1,2))
  
  Var_e1 =   T4* (cov(1,1)*Iyy*Iyy +cov(2,2)*Ixx*Ixx -2.0*Ixx*Iyy*cov(1,2))
  Var_e2 =   T4* (Ixy*Ixy*(cov(1,1)+cov(2,2)+2.0*cov(1,2)) -2.0*T*Ixy*(cov(1,3)+cov(2,3)) +T*T*cov(3,3))
  Cov_e1_e2= T4* (-Iyy*Ixy*cov(1,1)+Ixx*Ixy*cov(2,2)+(Ixx-Iyy)*Ixy*cov(1,2) + T*(Iyy*cov(1,3) -Ixx*cov(2,3))) 
  
  Corr=Cov_e1_e2/(sqrt(Var_e1*Var_e2))
  Var_T=Cov(1,1)+Cov(2,2)
  
  if keyword_set(print) then begin
      print
      print,'Var_e1',Var_e1
      print,'Var_e2',Var_e2
      print,'Cov_e1_e2',Cov_e1_e2
      print,'Corr',Corr
  endif
  
  return
end
