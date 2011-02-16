pro test_mom_error_es,numfr,aratio,tot
;test the moment errors

if n_params() eq 0 then begin
	print,'-syntax  test_mom_error_es,numfr,aratio,tot'
	return
endif

make_gaussian,psf,size=[10,10],fwhm=3.5,aratio=aratio,theta=0.0

arrsize = 2000
numobj = 50
spawn, 'rm /sdss/data1/esheldon/imcheck.fits'
for i=0, numfr-1 do begin
	im=fltarr(arrsize,arrsize)
	make_grid,arrsize,arrsize,50,50,x,y
;	im(x,y)=randomu(seed,2500)*55000.0+400.0
	im(x,y)=55000.0+400
	im=convol(im,psf,/edge_truncate)
	im=im+1000.0
	add_noise,im
	if i eq 0 then begin
		mwrfits,im,'/sdss/data1/esheldon/imcheck.fits'
	endif
;	extract,im,cat
	sdss_extract,im,cat
	if i eq 0 then begin
		tot=cat
	endif else begin
		concat_structs,tot,cat,temp
		tot=temp
	endelse
endfor

return
end
	

	

