pro sky,im,skyv,sigmav,silent=silent, niter=niter,nsig=nsig
;a replacement for the sky program which doesn't work very well
;just uses sigma clip

IF NOT keyword_set(silent) THEN silent = 0 ELSE silent =1
if n_params() eq 0 then begin	
	print,'sky,im,skyv,sigmav,silent=silent, niter=niter,nsig=nsig'
	print,'USING DAVES VERSION OF SKY'
	print,'USING DAVES VERSION OF SKY'	
	print,'USING DAVES VERSION OF SKY'

	return
end

IF n_elements(niter) EQ 0 THEN niter=2.0
IF n_elements(nsig) EQ 0 THEN nsig = 3.5
sigma_clip,im,skyv,sigmav,nsig=nsig,niter=niter, silent=silent

return
end


