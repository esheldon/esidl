function nfw_vel,r,p
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


Sqrt_G=65.583
;check this!

;del_c=(200/3.0)*c^3 /(alog(1.0+c)-c/(1.0+c))
v200=r200*Sqrt_G*sqrt(800.0*!pi*rho_crit/3.0)

T=alog(1.0+c*x)-c*x/(1.0+c*x)
B=alog(1+c)-c/(1.0+c)

return, v200*sqrt(T/(B*x))
end
