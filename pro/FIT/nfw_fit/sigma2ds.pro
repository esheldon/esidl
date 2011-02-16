;+
; NAME: sigma2ds
;
; PURPOSE:
; Computes DeltaSigma from 2D density Sigma
;
; CALLING SEQUENCE:
;    sigma2ds,r,sig,ds
;
; INPUTS:
; r  - array of binned radii
; Sig - the 2D density Sigma
;
; OUTPUTS:
;     ds - DeltaSigma
;
; METHOD:
;     Assumes power law interpolation as well as power law extrapolation
;-

PRO sigma2ds,r,sig,ds
IF n_params() EQ 0 THEN BEGIN
    print,'-syntax sigma2ds,r,sig,ds'
    return
endif

n=n_elements(r)
ds=dblarr(n)

al=alog(sig/shift(sig,1))/alog(r/shift(r,1))
al[0]=al[1]

slope_min=-1.95
if al[0] lt slope_min then begin
    print,'Warning - profile too steep to converge'
    print,'Truncating inner slope to ',slope_min
    al[0]=slope_min
endif

A=sig/(r^al)

RA=shift(r,1)
RA[0]=0.0
RB=r

Ints=A/(al+2.0) *(RB^(al+2.0)-RA^(al+2.0))
Icum=total(Ints,/cum)
avg_sig=2.0*Icum/(r^2)
ds=avg_sig-sig

return
end
