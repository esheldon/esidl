;+
; NAME:
;	icl_addobjects
;
; PURPOSE:
;	Add fake sersic profile(s) to an input image.
;
; CALLING SEQUENCE:
;	icl_addobjects, image, counts, sersicn, r50, xcen, ycen, 
;		pixscale=0.396, sblim=28.0, nx=, ny=
;
; INPUTS:
;	image: The original image
;	
;	The following describe the profile(s) to be added and must be all the
;	same length, either scalar or arrays.
;		counts: total counts in the image in nmgy.
;		sersicn: sersic index
;		r50:  Half-light radius in arcseconds.
;		xcen, ycen:  Location of center in the main input image.
;
; Optional Inputs
;	pixscale: default 0.396 arcsec.
;	sblim: Surface brightness limit in mag/square arcsec.
;	nx, ny: size of image.  NOTE the location of the center in the little
;		image is always nx/2 ny/2.  If these are not entered, the profile
;		is grown to he sblim surface brightness limit and added to the 
;		image.
;
; OUTPUTS:
;	The input image is modified.
;
; MODIFICATION HISTORY:
;	Created:  2008-11-05, Erin Sheldon, BNL
;
;-

pro icl_addobjects_syntax
	print,'icl_addobjects, image, counts, sersicn, r50, xcen, ycen, pixscale=0.396, sblim=28.0, nx=, ny='
	print,'  image: main image to which we will add'
	print,'  counts,... are arrays of sersic parameters.  xcen,ycen is the'
	print,'      center in the main image'
	message,'halting'
end

pro icl_addobjects, image, counts, sersicn, r50, xcen, ycen, $
		pixscale=pixscale, sblim=sblim, $
		nx=nx, ny=ny, slow=slow

	if n_params() lt 6 then begin
		icl_addobjects_syntax
	endif
	
	nadd=n_elements(counts)

	bad=0
	if n_elements(sersicn) ne nadd or $
			n_elements(r50) ne nadd or $
			n_elements(xcen) ne nadd or $
			n_elements(ycen) ne nadd then begin
		bad=1
	endif
	if n_elements(nx) ne 0 or n_elements(ny) ne 0 then begin
		if n_elements(nx) ne nadd or n_elements(ny) ne 0 then bad=1
	endif
	if bad then begin
		print,'All arrays must be the same length'	
		icl_addobjects_syntax
	endif

	if n_elements(pixscale) eq 0 then pixscale=0.396
	if n_elements(sblim) eq 0 then sblim=28.0

	sz=size(image, /dim)

	for i=0L, nadd-1 do begin

		if not keyword_set(slow) then begin

			; Calculate how big we must go to reach our surface brightness
			; limit for the input parameters
			if n_elements(nx) eq 0 then begin
				sbrad = sersic_sbrad(counts[i], r50[i], sersicn[i], sblim)
				; nx,ny twice the radius.   
				nxi = 2*sbrad/pixscale
				nyi = 2*sbrad/pixscale
				;print,'counts:',counts[i],'  r50:',r50[i],'  n:',$
				;	sersicn[i],'  sblim:',sblim
				;print,'Calculated size from surface brightness:',nxi
			endif else begin
				nxi = nx[i]
				nyi = ny[i]
			endelse

			im = dfakegal(				$
				flux=counts[i],			$
				sersicn=sersicn[i],		$
				r50=r50[i]/pixscale,	$ 
				nx=nxi,					$ 
				ny=nyi)

			;print,'flux in:'+ntostr(counts[i])+' flux out: '+ntostr(total(im))


			; location of columns and rows in big image coordinates

			; center in new image: nx/2L is how dfakegal finds the center
			imxcen = nxi/2L
			imycen = nyi/2L
			; the input xcen,ycen tells us where this goes in the main image.
			; so the 0,0 point of the image is this distance away
			col0 = xcen[i] - imxcen
			row0 = ycen[i] - imycen
			; now location of cols/rows in big image
			x = lindgen(nxi) + col0
			y = lindgen(nyi) + row0

			; find those pixels actually in the big image

			wx = where( x GT 0 AND x LT sz[0], nwx)
			wy = where( y GT 0 AND y LT sz[1], nwy)

			minwx=min(wx, max=maxwx)
			minwy=min(wy, max=maxwy)

			maxx=x[maxwx]
			minx=x[minwx]
			maxy=y[maxwy]
			miny=y[minwy]

			image[minx:maxx,miny:maxy] = $
				image[minx:maxx,miny:maxy] + im[minwx:maxwx, minwy:maxwy]
		endif else begin
			; simpler way
			im = dfakegal(				$
				flux=counts[i],			$
				sersicn=sersicn[i],		$
				r50=r50[i]/pixscale,	$ 
				xcen=xcen[i],			$
				ycen=ycen[i],			$
				nx=sz[0],				$ 
				ny=sz[1])

			image[*,*] = image[*,*] + im[*,*]

		endelse

	endfor
	
end
