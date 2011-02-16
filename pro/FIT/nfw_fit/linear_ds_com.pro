pro linear_ds_com,r,ds,_extra=ex,recache=recache,justcache=justcache

if n_params() eq 0 and not keyword_set(justcache)  then begin
    print,'-syntax linear_ds_com,r,ds,ns=ns,Omega_m=Omega_m,Omega_b=Omega_b,h=h,sigma_8=sigma_8'
    return
endif

;This assumes r is in comoving
;If not then call like
;IDL> linear_ds,r*(1+z),ds 
;IDL> growth_factor,1.0/(1+z),D,Omega_m
;IDL> ds=ds*D^2*(1.+z)^3

;call this once with an r that spans a long range
;example:
;IDL> rr=10.0^(dindgen(200)/199*4.3-2)
;IDL>linear_ds,rr,ds,/recache
;then for any r:
;IDL> linear_ds,r,ds
;so that it can interpolate from the cache

if keyword_set(justcache) then begin
    rr=10.0^(dindgen(200)/199*4.3-2)
    linear_ds_com,rr,junk,/recache
    return
endif


defsysv,'!linear_ds_cached',Ex=ex

if ex eq 1 and not keyword_set(recache) then begin
    ;print,'Using cached'
    rc=*!linear_r_cached
    dsc=*!linear_ds_cached
    ds=exp(interpol(alog(dsc),alog(rc),alog(r)))
    minr=min(rc,wmin)
    w=where(r lt minr,nums)
    if nums gt 0 then ds[w]=ds[wmin]
    return
endif

;num=200
;r=10.0^(dindgen(num)/(num-1)*3.1-1)

if ex then begin
    print,'freeing cached vars'
    ptr_free,!linear_xi_cached
    ptr_free,!linear_r_cached
    ptr_free,!linear_j3_cached
    ptr_free,!linear_sig_cached
    ptr_free,!linear_ds_cached
endif

;default params
if n_elements(ns) eq 0 then ns=0.95
if n_elements(Omega_m) eq 0 then Omega_m=0.27
if n_elements(Omega_b) eq 0 then Omega_b=0.045
if n_elements(h) eq 0 then h=0.71
if n_elements(sigma_8) eq 0 then sigma_8=1.0 

;num=200
;rr=10.0^(dindgen(num)/(num-1)*3.1-1.0)

linear_xi,double(r),xi,ns=ns,Omega_m=Omega_m,Omega_b=Omega_b,h=h,sigma_8=sigma_8,/pos

rho_crit=0.277520
rho_bar=rho_crit                ;assumes z=0

;keep in units of Omega_m*sigma_8^2
;do don't put in Omega_m, use rho_crit

drho=xi*rho_bar

drho2sigma,double(r),double(drho),sig
sigma2ds,double(r),double(sig),ds

xi2j3,r,xi,j3

print,'caching linear variables'
defsysv,'!linear_xi_cached',ptr_new(xi)
defsysv,'!linear_r_cached',ptr_new(r)
defsysv,'!linear_j3_cached',ptr_new(j3)
defsysv,'!linear_sig_cached',ptr_new(sig)
defsysv,'!linear_ds_cached',ptr_new(ds) 

return
end
