;+
; NAME:
;   sersic_sbrad
;
; PURPOSE:
;   Return the radius at which a 2-d sersic profile with the given input 
;   parameters drops to the specified surface brightness limit.
;
; CALLING SEQUENCE:
;   rlim = sersic_sbrad(flux_nmgy, r50, n, sblimit, /r0)
;
; INPUTS:
;   flux_nmgy: Flux in nanomaggies. Defined such that 
;		mag=22.5-2.5*alog10(flux_nmgy)
;   r50: The half-light radius.  Should be in arcsec if surface brightness 
;		limit is in mags/arcsec^2
;   n: Sersic index such that f=exp( -(r/r0)^(1/n))
;   sblimit: Surface brightness limit.  Should be in mags/arcsec^2 if r50 is in arcsec.
;
; KEYWORD PARAMETERS:
;   /r0: The input radius is r0 instead of r50.
;
; OUTPUTS:
;   radius at which the profile drops to the input flux density.
;		units will depend on units of r50 and surface brightness limit.
;
; MODIFICATION HISTORY:
;   Created: 2007-04-13, Erin Sheldon, NYU
;
;-
function sersic_sbrad, nmgy, r50_in, n, sblimit, r0=r0wasinput

    if n_params() lt 4 then begin
        print,'-Syntax: rad = sersic_sbradius(flux_nmgy, r50, n, sblimit, /r0)'
        on_error, 2
        message,'Halting'
    endif

    ; normalization of a 2-d sersic disk with these parameters
    if keyword_set(r0wasinput) then begin
        r0 = r50_in
        norm = sersic_norm(r0, n, /r0)
    endif else begin
        r0 = sersic_r0(r50_in, n)
        norm = sersic_norm(r0, n, /r0)
    endelse

    ; convert the requested surface brightness limit to a flux density limit
    fdlimit = icl_sb2fluxdensity(sblimit)

    ; solve for radius where this profile falls to that value
    arg = fdlimit*norm/float(nmgy)
    if arg gt 1.0 then begin
        rlim= -9999.0
    endif else begin
        rlim = r0*( -alog(arg) )^n
    endelse

    return, rlim 
end

