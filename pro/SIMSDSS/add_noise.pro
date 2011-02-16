pro add_noise,im,imout, gain=gain
;adds poisson noise (actually normally distributed noise) to an image
;if sky is less than 15 then it uses real poisson noise

if n_params() eq 0 then begin
	print,'-syntax  add_noise,im,imout, gain=gain'
        print,' gain only used if all values are larger that 15.0'
	return
endif

common seed,seed

s=size(im)
sx=s(1)
sy=s(2)

imout=fltarr(sx,sy)

small=where(im lt 15.0 AND im GT 0.0,nsmall)
big=where(im ge 15.0,nbig)
if ( not (keyword_set(gain) and nsmall eq 0) ) then gain=1.0

if nsmall gt 0 then begin 
print,'computing real poisson statistics for',strcompress(nsmall),' pixels'
	print,' -this is slow'
	imout(small)=poidev(im(small))
endif

if nbig gt 0 then begin
	imout(big)=im(big)+sqrt(1.0/gain*im(big))*randomn(seed,nbig)
endif

if n_params() eq 1 then im=imout

return
end
