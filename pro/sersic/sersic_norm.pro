;+
; NAME:
;   sersic_norm
;
; PURPOSE:
;   Calculate normalization of a 2-d sersic profile given the index n and r50 (or r0).
;   If maxr is set then the integral of the profile to that radius is calculated. If
;   r50 is entered it is converted to r0. The conversion is tested for n in [0.05,30] 
;   and may fail outside that range.
;
; CALLING SEQUENCE:
;   norm=sersic_norm(r50, n, aratio=, maxr=, /r0, status=)
;
; INPUTS:
;   r50: The half-light radius.
;   n: the index. 
; OPTIONAL INPUTS:
;   aratio: The axis ratio.
;   maxr: Max radius for integration in units of r50 unless /r0 in which case it
;       is in units of r0.
;   /r0: The input radius is r0 instead of r50.  No conversion is performed.
;
; OUTPUTS:
;   norm: The integral of exp( (r/r0)^(1/n) )*2*pi*r over the specified range.
; OPTIONAL OUTPUTS:
;   status: 0 for success, else for failure. Will only possibly fail if converting
;       r50 to r0.
; METHOD:
;   int( (r/r0)^(1/n) 2*pi*r ) = 2*pi*n*r0^2*aratio*gamma(2*n, maxt)
;       where maxt = (rmax/r0)^(1/n)
;       and gamma(2*n, maxt) = igamma(2*n,maxt)*gamma(2n)
;       This is because the IDL igamma(2*n,maxt) = gamma(2n,maxt)/gamma(2n)
;
; MODIFICATION HISTORY:
;   Created: 2007-04-11 Erin Sheldon, NYU 
;-


function sersic_norm, r50_in, n, aratio=aratio, maxr=maxr_in, r0=r0wasinput, status=status

    if n_elements(r50_in) eq 0 or n_elements(n) eq 0 then begin
        print,'-Syntax: norm = sersic_norm(r50, n, aratio=, maxr=, /r0, status=)'
        print,'If aratio present, r50 is assumed to be the largest axis'
        print,'maxr= is in units of r50'
        on_error, 2
        message,'Halting'
    endif

    ; deal with inputs
    if keyword_set(r0wasinput) then begin
        r0 = r50_in
        if n_elements(maxr_in) ne 0 then maxr=maxr_in
    endif else begin
        r0 = sersic_r0(r50_in, n, status=status)
        if status ne 0 then return, -1

        ; make sure maxr is in units of r0
        if n_elements(maxr_in) ne 0 then begin
            maxr = maxr_in*r0/r50_in
        endif
    endelse
    status = 0

    if n_elements(aratio) ne 0 then begin
        if aratio lt 0 or aratio gt 1 then begin
            on_error, 2
            message,'Aratio '+ntostr(aratio)+' out of range [0,1]'
        endif
    endif else begin
        aratio = 1.0
    endelse


    ; calculate normalization
    gam = gamma(2.0*n)
    if n_elements(maxr) ne 0 then begin
        ; igamma gives fraction of total to infinity
        maxt = maxr^(1.0/n)
        gam = igamma(2.0*n, maxt)*gam
    endif 

    return, 2.0*!pi*n*r0^2*aratio*gam
end


