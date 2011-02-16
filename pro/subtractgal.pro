pro  subtractgal, image, cat, imout, func=f

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;       
; PURPOSE:
;	
;
; CALLING SEQUENCE:
;      
;                 
;
; INPUTS: 
;
; INPUT KEYWORD PARAMETERS:
;
;       
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; CALLED ROUTINES:
; 
; PROCEDURE: 
;	
;	
;
; REVISION HISTORY:
;	
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if N_params() eq 0 then begin
     print,'-Syntax: '
     print,''
     print,'Use doc_library,"subtractgal"  for more help.'  
     return
  endif


imout = image
x0 = cat.x_image 
y0 = cat.y_image 

siz=size(image)
sx = long( siz[1] )
sy = long( siz[2] )

index=lindgen(sx*sy)
x = float(index mod sx)
y = float(index/sx)
r2 = sqrt( (x-x0)^2 + (y-y0)^2 )
maxvalue = 80.0

tol = .5
max=450

max= 500
rw = where(r2 le max-10)
xvals = indgen(max)
f = fltarr(max)
f[*] = image[x0,y0]
for i=1,max-1 do begin
   w = where(r2 ge i-1 and r2 le i+1)
   if (w[0] ne -1) then begin
     keep=[-1.0]
     for kk=0,n_elements(w)-1 do begin
       idx=index[ w[kk] ]
       if (image[idx] lt maxvalue) then begin
          if (keep[0] ne -1.0)  then keep=[keep,image[idx]] $
          else keep = [image[idx]]
       endif
     endfor
     f[i] =  f[i] + (keep[0] ne -1)*(median(keep) - f[i])
   endif
endfor

plot,f,psym=4
print,'Hit a key'
key=get_kbrd(20)   
if key eq 'q' then return


;   imout[index(rw)] = imout[index(rw)] - interpolate(f, r2[rw],cubic=-.5)

   imout[index(rw)] = imout[index(rw)] - interpol(f, xvals, r2[rw]) 
   imout[index(rw)] = imout[index(rw)] - $
           (imout[index(rw)] gt maxvalue)*imout[index(rw)]



rdis_setup,image,pls
rdis,image,pls

key=get_kbrd(20)

rdis_setup,imout,pls2
rdis,imout,pls2

return
end

   




























