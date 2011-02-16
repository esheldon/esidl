;+
; NAME: xi2sigma
;  
; PURPOSE: compute sigma(r) from xi(r)
;
; CALLING SEQUENCE: xi2sigma,r,xi,rr,sig
;
; INPUTS:
;      r - radii in Mpc/h
;      xi - 3D correlation function
;
; OUTPUTS:
;      rr - output radii in Mpc/h (not exactly the same as r) 
;      sig - sigma(rr)
;
; COMMENTS:
;
;  
; METHOD:
;      uses local power law interpolation
;      uses formula from Zehavi et al. 2005? 
;     (The intermediate scale clustering of LRGs)
; PROCEDURES CALLED:
; 
;
; INTERNAL SUPPORT ROUTINES 
;  
;
; REVISION HISTORY:
;   20-July-2005  Written by David Johnston, Princeton.
;  
;-
;------------------------------------------------------------------------------

pro xi2sigma,r,xi,rr,sig
;take linear xi versus r and make sigma
;uses local pl interpolation to perform integration
;uses Idit's formula
;output radius rr is just the r but truncated since one needs a factor
;of two larger in xi to compute sigma

if n_params() eq 0 then begin
    print,'-syntax xi2sigma,r,xi,rr,sig'
    return
endif

;!p.multi=[0,2,2]

wneg=where(xi lt 0,nneg)

if nneg gt 0 then begin
    print,'Error - xi must be positive for this powerlaw interpolation'
    return
endif

;choose a new set of points, logarithmically spaced such that a
;whole number of points correspond to a factor of two

nper2=100
;number of points per factor of 2 

n=n_elements(r)
rmax=max(r)
rmin=min(r)

num=1+nper2*ceil(alog(rmax/rmin)/alog(2))
rr=rmin*2.0^(dindgen(num)/nper2)
lrr=alog(rr)
lx=alog(xi)
lr=alog(r)
lxx=interpol(lx,lr,lrr)
xxi=exp(lxx)

;plot,r,xi,psym=1,/xlog,/ylog
;oplot,rr,xxi,psym=7,color=!magenta


al=(shift(lxx,-1)-lxx)/(shift(lrr,-1)-lrr)
al[num-1]=al[num-2]
wbad=where(al lt -2.95,nbad)

if nbad gt 0 then begin
    ;this doesn't need to be called , just can't be -3 exactly
    ;print,'Warning - log-slope of xi is sometimes < -3, something probably wrong'
    ;plot,rr,al,/xlog,ytit='al'
    ;return
endif

A=xxi/rr^(al)

;plot,rr,al,/xlog,ytit='al'
;plot,rr,A,/xlog,ytit='A'

Ra=rr
Ra[0]=0.0
;first one , integrate to zero
Rb=shift(rr,-1)

;3 integrals to do
beta=0.0
ex=3.0+al+beta
int0=A*(1.0/ex) *(Rb^ex -Ra^ex)

beta=1.0
ex=3.0+al+beta
int1=A*(1.0/ex) *(Rb^ex -Ra^ex)

beta=3.0
ex=3.0+al+beta
int2=A*(1.0/ex) *(Rb^ex -Ra^ex)

;array of sub integrals, last one to be ignored

;Now add up the sub integrals with 0 to rmin integral
T0=shift(total(int0,/cum),-nper2)
T1=shift(total(int1,/cum),-nper2)
T2=shift(total(int2,/cum),-nper2)

sig2=3.0*T0/(rr^3) + (-9.0/4)*T1/(rr^4) + (3.0/16)*T2/(rr^6)

;now trim off last bunch
rr=rr[0:num-nper2-2]
lrr=lrr[0:num-nper2-2]
sig2=sig2[0:num-nper2-2]
sig=sqrt(sig2)

;plot,rr,sig,/xlog,/ylog,ytit='sigma',psym=1,yr=[min(sig)*0.9,max(sig)*1.1],/yst
;wtest=where(rr gt .01 and rr lt .3 and (indgen(num) mod fix(nper2/2)) eq 0)
;rrr=rr[wtest]
;xi2sigma_qromb,rr,xxi,rrr,sss
;oplot,rrr,sss,color=!magenta,psym=7
;plot,rr[wtest],sig[wtest]/sss,/xlog,/yno,psym=-1

;now interpolate back onto original points

ss=exp(interpol(alog(sig),lrr,lr))
wcut=where(r*2 lt rmax) 
rr=r[wcut]
sig=ss
sig=sig[wcut]

return
end


