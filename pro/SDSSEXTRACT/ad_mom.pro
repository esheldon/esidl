pro ad_mom,cat, im, ixx, iyy, ixy, err, numiter, wcenx, wceny,$
           whyflag,rho4,sky=sky,dontadd=dontadd

;this is the IDL wrapper for calling ad_mom.f
;ad_mom.f is a fortran program written by Phil Fischer
;which calculates Adaptively Weighted Moments
;it adjusts the weighting function until the measurement 
;is optimal. This happens when the moments of the weight 
;are twice the weighted moments. Then the unweighted moments 
;are just equal to the moments of the weight and therefore 
;twice the weighted moments. It returns the unweighted moments 
;measured in this weighted way.
;cat: the sextractor catalog (an IDL structure)
;im: the image
;ixx,iyy,ixy the output adaptive moments
;err: a measure of uncertainty is 9.99 if moments are
;not found (ie. does not converge sensibly)
;numiter: the maximum number of iterations
;wcenx,wceny: the output weighted centroids
;whyflag: the reason why it did not converge 
;see the fortran code for their meanings
;sky:input an array of sky values for each objects
;if not present it will calculate the sky for the whole image
;and use this one number for all  
;if im is undefined it will run on a test image
;dontadd: if this keyword is set it will NOT add one to
;the x and y coordinates before handing to fortran
;the default is to assume that the x,y are already in
;IDL notation that is the first pixel in (0,0)
;If however the catalog comes directly out of sextractor
;it will already be in the same notation as fortran (ie. 
;first pixel is (1,1)) so set the keyword dontadd
; -Dave Johnston
; -Erin Scott Sheldon  Changed sofile, added comments
;       changed sky stuff

if n_params() eq 0 then begin
	print,'-syntax ad_mom,cat,im,ixx,iyy,ixy,err,'
	print,'numiter,wcenx,wceny,rho4,'
	print,'whyflag,sky=sky,dontadd=dontadd'
	return
endif

im=float(im)

ixx=float(cat.x2_image)
iyy=float(cat.y2_image)
ixy=float(cat.xy_image)

x=float(cat.x_image)
y=float(cat.y_image)

if keyword_set(dontadd) eq 0 then begin
	x=x+1.0
	y=y+1.0
endif 

mag=float(cat.mag_best)

n=n_elements(x)
sh=2.0

IF n_elements(sky) EQ 0 THEN BEGIN 
  sky,im,sk,smode, sigsk
  sigsk=replicate(sigsk,n)
  sk=replicate(sk,n)
ENDIF ELSE BEGIN
  sk=sky
  sigsk=5.9
  sigsk=replicate(sigsk,n)
ENDELSE

numiter=lonarr(n)
wcenx=fltarr(n)
wceny=fltarr(n)
whyflag=lonarr(n)

s=size(im)
nx=s(1)
ny=s(2)
err=x
rho4 = fltarr(n)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; define the shared object file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

bindir='/sdss3/usrdevel/philf/idl.lib/'
sofile = bindir + 'ad_momi.so'
keyword = 'ad_mom_'

print
print,'Using sofile:  ',sofile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Call the fortran routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print,'Calling ad_momi.f'
ff=call_external(sofile, keyword, x, y, ixx, iyy, ixy, n, sh, im, sk, nx, ny,$
                 err, sigsk, mag, numiter, wcenx, wceny, whyflag, rho4)


return
end























