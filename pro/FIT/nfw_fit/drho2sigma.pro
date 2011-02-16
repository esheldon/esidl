;+
; NAME: drho2sigma
;  
; PURPOSE:
; Computes 2D density Sigma from delta_rho
;
; CALLING SEQUENCE:
;    drho2sigma,r,drho,sig,/extrap
;
; INPUTS:
; r  - array of binned radii, supposed to be the radius such that
;      drho=DeltaRho(r)  
; drho - DeltaRho
;  
; OPTIONAL KEYWORD INPUTS:
;      extrap - performs an endpoint correction (deafault)
;
; OUTPUTS:
;     sigma - the 2D density
;  
; METHOD:
;     Uses three-point neighborhoods to interpolate the function as a
;     quadratic and uses analytic formulas to compute the integral
;     endpoint correction computed by extrapolating the function and
;     recursively calling it self
;
; PROCEDURES CALLED:
;  laq_quad_poly - compute the Lagrangian interpolation quadratic polynomial
;-


pro drho2sigma,r,drho,sig,extrap=extrap
;Do the projection with tabulated data
;include extrap and it will entend the data w PL extrapolation to
;do an endpoint correction

if n_params() eq 0 then begin
   print,'-syntax  drho2sigma,r,drho,sig,extrap=extrap '
   return
endif

n=n_elements(r)

if n_elements(extrap) eq 0 then extrap=1 ;the default is to correct

If keyword_set(extrap) then begin
    fac=1.0                     ;the fraction of the log interval to add on, 0.5 is probably good
    next=long(n*fac > 10L)
    Lext=(alog10(max(r))-alog10(min(r)))*fac
    r_ext=max(r)*10.0^(Lext*(dindgen(next)+1)/(next-1))
    al=alog(drho[n-1]/drho[n-2])/alog(r[n-1]/r[n-2])
    A=drho[n-1]/(r[n-1]^al)
    drho_ext=A*r_ext^al
    
    r_ext=[r,r_ext]
    drho_ext=[drho,drho_ext]
    drho2sigma,r_ext,drho_ext,sig,extrap=0
    Sig=sig[0:n-1]
    return
endif

Sig=dblarr(n)

for j=0, n -3 do begin

    RR=r[j]
    num=n-j
    w=lindgen(num)+j

    Int=dblarr(num)
    
    for i=0L, num-3 do begin
        x=r[i+j:i+j+2]
        y=drho[i+j:i+j+2]
        p=lag_quad_poly(x,y)
        A=x[0]
        B=x[1]
        
        Rad=B    
        S=sqrt(Rad^2-RR^2)
        I0_B = S
        I1_B = S*Rad/2.0 + RR^2*alog(Rad+S)/2.0
        I2_B = S*(2*RR^2 + Rad^2)/3.0
        
        Rad=A    
        S=sqrt(Rad^2-RR^2)
        I0_A = S
        I1_A = S*Rad/2.0 + RR^2*alog(Rad+S)/2.0
        I2_A = S*(2*RR^2 + Rad^2)/3.0
    
        I0=(I0_B-I0_A)*p[0]
        I1=(I1_B-I1_A)*p[1]
        I2=(I2_B-I2_A)*p[2]
        
        Int[i]=2*(I0+I1+I2)
    endfor
    
    sig[j]=total(Int)
endfor

return
end

