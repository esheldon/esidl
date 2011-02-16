;+
; NAME:
;   tvellipse_mom
;
; PURPOSE:
;   Plot ellipses on objects based on their second moments, image must
;   be displayed first
;
; CALLING SEQUENCE:
;   tvellipse_ellipse,x,y,Ixx,Iyy,Ixy,scale=scale,_extra=ex
;
; INPUTS:
;   x,y -  centroids , can be arrays
;   Ixx,Iyy,Ixy - second moments (usually weighted) of the light about the centroids 
;
; OPTIONAL KEYWORDS
;   scale - scales the overall size of the ellipses by some number
;           default 1. (can be array)
;
; OUTPUTS: none
;
; COMMENTS:
;   This procedure checks to make sure the moments are sane. Only sane
;   objects gets ellipses. Insane objects get red circles. For example
;   Ixx,Iyy must be > 0 and and Ixy must be < sqrt(Ixx*Iyy).
;
; EXAMPLES:
;
; BUGS: 
;
; PROCEDURES CALLED:
;   tvellipse
;
; REVISION HISTORY:
;   09-May-2005  Written by David Johnston,Princeton
;-
;------------------------------------------------------------------------------

pro tvellipse_mom,x,y,Ixx,Iyy,Ixy,scale=scale,_extra=ex

if n_params() eq 0 then begin
    print,'-syntax tvellipse_ellipse,x,y,Ixx,Iyy,Ixy,scale=scale'
    return
endif

if n_elements(scale) eq 0 then scale=1.0
n=n_elements(x)

;do sanity checks
Imin=0.01
Tmax=2500.0
wgood=where(Ixx gt Imin and Iyy gt Imin and Ixy^2 lt Ixx*Iyy $
and Ixx+Iyy lt Tmax,ngood)
good=bytarr(n)
good[wgood]=1

T=(Ixx+Iyy) > 2*Imin

e1=(Ixx-Iyy)/T
e2=2*Ixy/T
e=sqrt((e1^2+e2^2)>0.0)
r=sqrt(((1-e)/(1+e))>0.0)

hard_scale=4.0                  ;scaling their sizes?
S=sqrt(T)*hard_scale*scale
Rmax=S
Rmin=S*r
theta=(180.0/!pi)*0.5*atan(2*Ixy,(Ixx-Iyy))

for i=0, n-1 do begin
    if good[i] then begin
        ;good ones in yellow
        tvellipse,Rmax[i],Rmin[i],x[i],y[i],Theta[i],color=!yellow,/data,_extra=ex
    endif else begin
        ;bad ones in red
        tvcircle,10,x[i],y[i],/data,color=!red
    endelse
endfor

return
end
