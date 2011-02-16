;+
; NAME:
;   sersic_image
;
; PURPOSE:
;   Returns a normalized 2-d sersic profile projected onto a pixel array.  The
;   normalization is to infinity, so the total flux in the image may be less. 
;	The profile can be elliptical, oriented at any angle, and off-center.
;
; CALLING SEQUENCE:
;   im = sersic_image(r0, n, [xsize,ysize], counts=1, aratio=1, 
;                       theta=0, cen=imagecen, /r50, status=)
;
; INPUTS:
;   r50: The sersic half-light radius.  Will be converted to r0.  This 
;		conversion has been tested to work for n ranging [0.05,20].
;   n: The Sersic index.
;   [xsize,ysize]: An array containing the image size.
;
; OPTIONAL INPUTS:
;   counts: The counts in the profile if it were integrated to infinity.
;   aratio: axis ratio
;   theta: orientation.
;   cen: Image center.  Default is (nx-1)/2, (ny-1)/2
;   /r0: The input radius is r0 instead of r50 and must be converted to r0. 
;
; OUTPUTS:
;   image
; OPTIONAL OUTPUTS:
;   status: 0 for success, else for failure. Will only possibly fail if 
;		converting r50 to r0.
;
; MODIFICATION HISTORY:
;   Created: 2007-04-11 Erin Sheldon NYU
;
;-

function sersic_image, r50_in, n, siz, counts=counts, aratio=aratio,$
                 theta=theta, cen=cen, r0=r0wasinput, core=core, status=status


    if n_params() LT 3 then begin
        print, '-syntax: image=sersic_image(r50, n, [xsize,ysize], counts=1, aratio=1, theta=0, cen=, /r0, status=)'
        print,'  Returns normalized 2-d sersic profile. Theta in degrees'
        on_error,2
        message,'Halting'
    endif

    if keyword_set(r0wasinput) then begin
        r0 = r50_in
        status=0
    endif else begin
        r0 = sersic_r0(r50_in, n, status=status)
        if status ne 0 then return, -1
    endelse
 
    ; Note: norm is only the real norm for aratio=0.0
    norm = sersic_norm(r0, n, /r0, aratio=aratio)

    sx=long(siz[0])
    sy=long(siz[1])

    ;; check some parameters

    if n_elements(counts) eq 0 then counts = 1.0 else counts=float(counts)
    if n_elements(aratio) eq 0 then aratio = 1.0 else aratio = float(aratio)
    if n_elements(theta) eq 0 then theta = 0.0 else theta = float(theta)
    if n_elements(cen) eq 0 then cen = [(sx-1.0)/2., (sy-1.)/2.]

    cx = float(cen(0))
    cy = float(cen(1))

    ;; create our galaxy
    image=fltarr(sx,sy)

    index= lindgen(sx*sy)
    x=index mod sx
    y=index/sx

    ct=cos(theta*!Pi/180.)
    st=sin(theta*!Pi/180.)

    ;; Rotate the coordinate system
    xp=(x-cx)*ct + (y-cy)*st
    yp=(-1*(x-cx)*st + (y-cy)*ct)/aratio

    ;; define r in this coordinate system
    r=sqrt(xp^2+yp^2)

    ;; scale r
    if n_elements(core) ne 0 then begin
        arg = ( (r+core)/r0 )^2
    endif else begin
        arg = (r/r0)^2
    endelse
    arg = arg^(1.0/(2*n))

    w=where(arg LT 87.0, good)    ;watch floating point underflow
    IF good GT 0 THEN BEGIN
        image[index[w]] = exp(-arg[w])/norm*counts
    ENDIF

    return, image
end


