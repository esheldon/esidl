function nfw_sigma,r,p
;p[0] is r200 the virial radius
;p[1] is c the concentration parameter

defsysv,'!nfw_200',ex=ex
if ex eq 0 then begin
    print,'ERROR - Must define !nwf_200 with define_nfw_200'
    return,-1
endif

r200=p[0]
c=p[1]
x=r/r200

rho_crit=0.277520
;Don't assumes at z=0,check for glbal vars
z=!nfw_200.z
Omega_m=!nfw_200.Omega_m
;z dependence of rho_crit from Freidman Eq, assuming flat LCDM
H2=Omega_m *(1.0 +z)^3 + (1-Omega_m)
rho_crit=rho_crit*H2

del_c=(200/3.0)*c^3 /(alog(1.0+c)-c/(1.0+c))

;rho=rho_crit*del_c /(c*x * (1.0+c*x)^2)
xx=double(c*x)
x=0.0
;redefined to the Wright and Brainerd
rs=r200/c
fac=2*rs*del_c*rho_crit

;three cases to consider
ep=double(0.001)                          ;buffer , use linear interp in here

sig=dblarr(n_elements(r))
w1=where(xx lt (1.0-ep),n1)
w2=where(xx ge (1.0-ep) and xx le (1.0+ep),n2)
w3=where(xx gt (1.0+ep),n3)

if n1 gt 0 then begin
    x=xx[w1]
    s=(1 -2/sqrt(1 - x^2)*atanh(sqrt((1-x)/(1+x))))/(x^2-1)
    sig[w1]=s
endif

if n3 gt 0 then begin
    x=xx[w3]
    s=(1 -2/sqrt(x^2 -1)*atan(sqrt((x-1)/(1+x))))/(x^2-1)
    sig[w3]=s
endif

if n2 gt 0 then begin
    x=xx[w2]
    s=1/3d                     
    x1=1.0-ep
    x2=1.0+ep
    s1=(1 -2/sqrt(1 - x1^2)*atanh(sqrt((1-x1)/(1+x1))))/(x1^2-1)
    s2=(1 -2/sqrt(x2^2 -1)*atan(sqrt((x2-1)/(1+x2))))/(x2^2-1)
    s=(x-x1)*s2/(x2-x1) + (x-x2)*s1/(x1-x2)
    sig[w2]=s
endif

sig=sig*fac

return,sig
end
