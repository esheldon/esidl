;+
; NAME:   linear_xi
;  
; PURPOSE: Compute the linear correlation function using Eisenstein
;          and Hu Transfer function (can be generalized) 
; 
; CALLING SEQUENCE:
;          linear_xi,r,xi,ns=ns,Omega_m=Omega_m,Omega_b=Omega_b,h=h,sigma_8=sigma_8,$
;           bbks=bbks,plot=plot,_extra=ex
;
; INPUTS:
;        r - radius is Mpc/h, can be array
; 
; REQUIRED KEYWORD INPUTS:
;        ns - scalar spectral index
;        Omega_m - matter density parameter, including baryons
;        Omega_b - baryon density parameter
;        h - hubble constant in units of 100 km/s/Mpc
;        sigma_8 - normalization 
;
; OPTIONAL KEYWORD INPUTS:
;        bbks - if set will use BBKS transfer instead of Eisenstein and Hu
;        plot - if set will make a plot
;        _extra - extra plotting keywords
;
; OUTPUTS:
;        xi - the 3D linear correlation function unitless
;
; COMMENTS:
;  
; METHOD:
;         does local powerlaw interpolations for integrals when Fourier transforming
;
; PROCEDURES CALLED:
;         eisen_hu - computes transfer function
;         transfer_bbks - optionally use this transfer
;         xi_integral - function for performing the integrations
;         xi2sigma - needed for normalization
;
; REVISION HISTORY:
;   20-July-2005  Written by David Johnston, Princeton.
;  
;-
;------------------------------------------------------------------------------

pro linear_xi,r,xi,ns=ns,Omega_m=Omega_m,Omega_b=Omega_b,h=h,sigma_8=sigma_8,bbks=bbks,$
    nplk=nplk, $
    cmbfast=cmbfast,$
    tfile=tfile,plot=plot,$
    _extra=ex,pos=pos,$
    k=k,Pk=Pk,nonorm_pk=nonorm_pk
;calculate the linear correlation function
;uses Eisenstein and Hu unless bbks is set

if n_params() eq 0 then begin
    print,'-syntax linear_xi,r,xi,ns=ns,Omega_m=Omega_m,Omega_b=Omega_b,h=h,sigma_8=sigma_8'
    return
endif

;default params
if n_elements(ns) eq 0 then ns=0.98
if n_elements(Omega_m) eq 0 then Omega_m=0.273
if n_elements(Omega_b) eq 0 then Omega_b=0.0546
if n_elements(h) eq 0 then h=0.663
if n_elements(sigma_8) eq 0 then sigma_8=1.0 
if n_elements(tfile) eq 0 then tfile='~/data/cmbfastTF/SanchezTF.fits'

if n_elements(nplk) eq 0 then nplk=20.0

n=n_elements(r)
lkmin=-5.0d
lkmax=6.0d
numk=(lkmax-lkmin)*nplk
k=10.0^(dindgen(numk)/(numk-1.0)*(lkmax-lkmin)+lkmin)

fakenorm=5d6
rmax=2000.0

if max(r) gt rmax then begin
    print,'Warning some r larger than rmax, rmax=',rmax,'max(r)=',max(r) 
endif

;k=10.0^(dindgen(1000)/90-3.5)
if keyword_set(cmbfast) then begin
    print,'using CMBfast: ',tfile
    cmb=mrdfits(tfile,1)
    ns=cmb.n
    Tcmb=cmb.T
    Tcmb=Tcmb/Tcmb[0]
    kcmb=cmb.k
    kcmbmax=max(kcmb)
    T=exp(interpol(alog(Tcmb),alog(kcmb),alog(k))) 
    uu=((k-kcmbmax)/(5.0*kcmbmax))*(k gt kcmbmax) < 5.0
    Filt=exp(-uu)
    ;truncate it at high k?
    ;T=T*Filt
endif else begin
    if keyword_set(bbks) then begin
        transfer_bbks,k,T,Omega_m,h,Omega_b
    endif else begin
        eisen_hu,k,T,Omega_m=Omega_m,Omega_b=Omega_b,h=h
    endelse
endelse

Pk=fakenorm*k^ns * T^2

xi=dblarr(n)
for i=0L, n-1 do begin
    rr=r[i] < rmax
    ;xi_integral,k,Pk,rr,xii,d,plot=plot,yr=[.9999,1.0001],_extra=ex
    xi_int_pk,k,Pk,rr,xii,doplot=doplot
    if finite(xii) eq 0 then begin
        print,'Error-Nan   r=',rr
        stop
    endif
    xi[i]=xii
    if keyword_set(plot) then begin
        print,'OK?'
        key=get_kbrd(20)
        if key eq 'q' then return
        if key eq 's' then stop
    endif
endfor

;print,'Fake Normalization'
;make this faster

w=where(r lt 25.0)
xi2sigma,r[w],xi[w],rr,ss

s8=exp(interpol(alog(ss),alog(rr),alog(8.0)))

xi=xi*(sigma_8/s8)^2
nonorm_pk=Pk
Pk=Pk*(sigma_8/s8)^2

if keyword_set(pos) then begin
    ;make it positive
    ;make it drop off like a steep powerlaw
    ;useful for programs which will extrapolate it with a powerlaw
    beta=6.0
    dfac=0.3
    w=where(xi lt 0,nneg)
    if nneg gt 0 then begin
        w0=w[0]
        print,'Some negative ',nneg
        rr=r[w0-3:*]
        xx=xi[w0-3:*]
        ;plot,rr,xx,psym=-1
        ;oplot,rr,rr*0,color=!red,linest=2
        A=xi[w0-1]
        RA=r[w0-1]
        xx2=dfac*A*(rr/RA)^(-beta)
        xi[w]=dfac*A*(r[w]/RA)^(-beta)
        ;oplot,rr,xx2,color=!magenta,psym=-7
    endif
endif

return
end
