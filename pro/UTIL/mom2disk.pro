;+
; NAME:
;	mom2disk
;
; PURPOSE:
;	Take input Ixx,Ixy,Iyy and return an image with pixel values corresponding
;	to those values.  The model can be 'gauss','exp' or 'dev'
;
;       Note the moments don't really translate properly for exp
;       and for dev, but this is just a way to use the same sorts of inputs.
;
;
; CALLING SEQUENCE:
;	gauss=mom2disk(model, Ixx, Ixy, Iyy, imsize, counts=, cen=, status=)
;
; INPUTS:
;   model:'gauss' 'exp' or 'dev'
;	Ixx, Ixy, Iyy: Moments of the distribution. These could be the moments
;		returned by adaptive moments for example.
;	imsize: Size of the image in the x and y directions [sx, sy]
;
; OPTIONAL INPUTS:
;	counts:  Default is 1.0
;	cen: Default is (sx-1)/2. and (sy-1)/2.
;
; OUTPUTS:
;	An image of size [sx,sy].
;
; MODIFICATION HISTORY:
;	Created 2010-02-13, Erin Sheldon, BNL
;
;-
function mom2disk, model, ixx, ixy, iyy, imsize, $
               counts=counts, cen=cen, status=status

    status = 1
    if n_params() LT 5 then begin 
        print,'-Syntax gauss = mom2disk(model, ixx, ixy, iyy, imsize, counts=, cen=)'
        print,'produce an image from the input ixx,ixy,iyy'
        on_error, 2
        message,'halting'
    endif

    det = ixx*iyy - ixy^2
    if det eq 0.0 then begin
        message,'Error:  determinant is zero'
    endif

    sx=long(imsize[0])
    sy=long(imsize[1])

    if n_elements(counts) eq 0 then counts=1.
    if n_elements(cen) eq 0 then begin
        cx = (sx-1.)/2.
        cy = (sy-1.)/2.
        cen=[cx,cy]
    endif else begin
        cx=cen[0]
        cy=cen[1]
    endelse


    ; Compute the exponent of the gaussian at each point in the image.  
    ; Use clever indexing technique

    index=lindgen(sx*sy)

    x=index MOD sx
    y=index/sx

    rr = (x-cx)^2*Iyy -2*(x-cx)*(y-cy)*Ixy + (y-cy)^2*Ixx

    case strlowcase(model) of
        'gauss': begin
            rr = 0.5*rr/det
        end
        'exp': begin
            rr = sqrt(rr/det)
        end
        'dev': begin
            rr = sqrt(rr)
            rr = 7.67*( (rr/det)^(.25) -1 )
        end
        else: message,'model must be one of gauss, exp, or dev'
    endcase


    ; Compute the gaussian
    disk=fltarr(sx,sy)
    w=where(rr lt 10.8, nw)      ; watch floating point underflow
    ; this will allow ~.00002 as smallest number
    ; in disk

    if nw gt 0 then disk[index[w]]=exp(-rr[w])

    ;; set the counts
    norm=total(disk)
    if norm gt 0 then begin
        disk = counts/total(disk)*disk
    endif
    return, disk


end
