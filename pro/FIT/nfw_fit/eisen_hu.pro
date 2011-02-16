;+
; NAME: eisen_hu
;  
; PURPOSE: compute the Eisenstein and Hu Transfer function
; 
;
; CALLING SEQUENCE: eisen_hu,k,T,Omega_m=Omega_m,Omega_b=Omega_b,h=h,plot=plot
;
;
; INPUTS: k - wavenumber in h/Mpc
;   
; REQUIRED KEYWORD INPUTS:
;         Omega_m - matter density parameter, including baryons
;         Omega_b - baryon density parameter
;         h - hubble constant in units of 100 km/s/Mpc
; 
; OPTIONAL KEYWORD INPUTS:
;         plot - will make a plot of the transfer function
;
; OUTPUTS:
;         T - the matter transfer function
;      
; COMMENTS:
;         tested against Bullocks code, they agree to better than 1% 
;         The Temp of the CMB is not specified in Eisenstein and Hu,
;         does it matter?
; METHOD:
;         Applies the formulas of Eisenstein and Hu 1998 , ApJ 496,605
;
; REVISION HISTORY:
;   23-July-2005  Written by David Johnston, Princeton.
;  
;-
;------------------------------------------------------------------------------

pro eisen_hu,kk,T,Tc,Tb,Omega_m=Omega_m,Omega_b=Omega_b,h=h,plot=plot,stil=stil
;Eisenstein and Hu transfer function
;Equations 2-7, 10-12, 14-24
;kk is in units of h/Mpc

if n_params() eq 0 then begin
    print,'-syntax eisen_hu,kk,T,Tc,Tb,Omega_m=Omega_m,Omega_b=Omega_b,h=h,plot=plot,ns=ns,Pk=Pk'
    return
endif

k=kk*h                          ;now in units 1/Mpc physical wavenumbers

Tcmb=2.728d                      ;look up more recent
Tcmb=2.72511d   ;What Bullocks code uses 

Theta=Tcmb/2.7d
frac=double(Omega_b)/Omega_m            ;baryon fraction     

Om=double(Omega_m)*h^2
Ob=double(Omega_b)*h^2

;done with h dependence

;Equation 2
Zeq=2.50d4 * Om * Theta^(-4)

;Equation 3
Keq=7.46d-2 * Om * Theta^(-2)

;Equations 4
b1=0.313d * Om^(-0.419) * (1d + 0.607 * Om^0.674)
b2=0.238d * Om^0.223d
Zd= 1291d * Om^0.251d * (1d + b1 * Ob^b2) / (1+0.659d * Om^0.828d)

;Equation 5

Rconst=31.5d * Ob * Theta^(-4)
Rd  = Rconst * (Zd/1e3)^(-1)
Req = Rconst * (Zeq/1e3)^(-1)

;Equation 6

s=2d/(3d*Keq) * sqrt(6d/Req) * alog((sqrt(1+Rd)+sqrt(Rd+Req))/(1+sqrt(Req)))
;Equation 7

Ksilk = 1.6d * Ob^0.52 * Om^0.73 * (1d + (10.4d * Om)^(-0.95))

;Equations 11,12

a1 = (46.9d * Om)^0.670 * (1d + (32.1 * Om)^(-0.532))
a2 = (12.0d * Om)^0.424 * (1d + (45.0 * Om)^(-0.582))
alc=a1^(-frac) * a2^(-frac^3) 

bb1 = 0.944d * (1d + (458.0d * Om)^(-0.708))^(-1)
bb2 = (0.395d * Om)^(-0.0266)
betac = (1d + bb1*((1.0d -frac)^bb2-1d))^(-1)

;Equation 10

q=k/(13.41d * Keq)

;Equations 18-20

e=exp(1d)
nat = alog(e + 1.8 * betac * q)
Term1=386d/(1d + 69.9d * q^1.08)
Co= 14.2/alc + Term1
To  = nat/(nat+ Co * q^2)
f=1.0/(1.0+(k*s/5.4)^4)
Co1 = 14.2d + Term1
To1 = nat/(nat+ Co1 * q^2)

;Equation 17

Tc= f * To1 + (1.0-f) * To

;Now the baryons

;Equation 14,15


y=(1.0d + Zeq)/(1d + Zd)
sq=sqrt(1d + y)
Gy=y * (-6d*sq + (2+3*y)*alog((sq+1.0)/(sq-1.0)))
alb=2.07d*Keq*s*(1d + Rd)^(-0.75)*Gy

;Equations 22,23

betanode = 8.41 * Om^0.435
stil=s*(1d + (betanode/(k*s))^3)^(-0.3333333)

;Equation 24

betab=0.5d + frac + (3d -2.0d *frac) * sqrt(1d +(17.2*Om)^2)
nat2 = alog(e+1.8*q) 
To11=nat2/(nat2+Co1*q^2)

x=k*stil
;klam=2*!pi/stil
;print,'klam ',klam
;stop
jo=sin(x)/x
Tb=(To11/(1d +(k*s/5.2)^2) + alb/(1.0+(betab/(k*s))^3)*exp(-(k/Ksilk)^1.4))*jo

T=frac*Tb + (1d -frac)*Tc

if keyword_set(plot) then begin
    plot,k,T^2,psym=3,/xlog,/ylog
    oplot,k,(Tb*frac)^2,psym=3,color=!blue
    oplot,k,(Tc*(1.0-frac))^2,color=!red
endif

return
end
