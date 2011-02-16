pro test_mom_error,numfr,tot
;test the moment errors

if n_params() eq 0 then begin
	print,'-syntax  test_pt,numfr,tot'
	return
endif

make_gaussian,psf,size=[7,7],fwhm=2.5,aratio=.8,theta=0.0

for i=0, numfr-1 do begin
	im=fltarr(2000,2000)
	make_grid,2000,2000,50,50,x,y
	im(x,y)=randomu(seed,2500)*55000.0+400.0
	im=convol(im,psf,/edge_truncate)
	im=im+1000.0
	add_noise,im
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
	

	

