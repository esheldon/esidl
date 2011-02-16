function admom::init

	common _tmblock, tmconv, tm_in_regauss
	tmconv=0.0
	tm_in_regauss=0.0

	; we should cut harder on this later I think
	self.detf0_lim = 0.1
	; this is absolute limit
	self.detf0_tol = 1.e-5
	return, 1
end

function admom::detf0_lim
	return, self.detf0_lim
end


;
; Run regaussianization on a set of objects from "str"
;
; To run a single object, see ::do_regauss
;

function admom::regauss, str, $
		indices=indices, $
		bands=bands, $
		unload=unload, $
		verbose=verbose

	; str must have these tags:
	;   run,rerun,camcol,field,id
	;   colc,rowc
	; If the tag m_rr_cc exists it will be used for a first guess at the
	; object size

	common admom_psfield_block_full, run, rerun, camcol, field, psfield, psp
	common _tmblock, tmconv, tm_in_regauss

	tm0=systime(1)
	tmread=0.0
	tmother=0.0
	tmregauss=0.0
	tmconv=0.0
	tm_in_regauss=0.0

	; we will return a struct for all, but may only process a subset
	; or none through admom
	ntotal = n_elements(str)
	if n_elements(indices) eq 0 then begin
		nprocess=ntotal
		indices=lindgen(nprocess)
	endif else begin
		if indices[0] eq -1 then nprocess=0 else nprocess=n_elements(indices)
	endelse

	if n_elements(bands) eq 0 then bands=[0,1,2,3,4]
	nband=n_elements(bands)

	outstruct = self->regauss_struct(ntotal)
	outstruct.run=str.run
	outstruct.rerun=str.rerun
	outstruct.camcol=str.camcol
	outstruct.field=str.field
	outstruct.id=str.id

	for i=0L, nprocess-1 do begin
		index=indices[i]

		if keyword_set(verbose) then begin
			i1=index +1
			if (nprocess gt 1) and (i1 mod 100) eq 0 then begin
				print,'Processing ',index+1,nprocess,$
					format='(a,i0,"/",i0)'
			endif
			if verbose gt 1 then begin
				print,$
					str[index].run, $
					str[index].rerun,$
					str[index].camcol, $
					str[index].field, $
					str[index].id
			endif
		endif

		tm1=systime(1)
		read_atlas, $
			str[index].run, $
			str[index].camcol, $
			str[index].field, $
			str[index].id, $
			rerun=str[index].rerun, $
			atlas_struct=astr, $
			status=status, /silent
		tmread += systime(1)-tm1

		if status eq 0 then begin

			tm1=systime(1)
			self->memcache_psfield_full, str, index
			tmread += systime(1)-tm1

			for iband = 0, nband-1 do begin

				band=bands[iband]

				psf = sdss_psfrec($
					*psp[band], $
					str[index].rowc[band], str[index].colc[band],$
					counts=1.0, /trim)

				tm1=systime(1)
				self->do_regauss_atlas, $
					str, outstruct, index, band, astr, psf
				tmregauss += systime(1)-tm1

			endfor ; loop over bandpasses
		endif ; success reading atlas
	endfor ; loop over objects

	print,'total ',f='(a,$)'
	ptime,systime(1)-tm0
	print,'read: ',f='(a,$)'
	ptime,tmread
	print,'regauss: ',f='(a,$)'
	ptime,tmregauss
	print,'convol: ',f='(a,$)'
	ptime,tmconv
	print,'inside internal regauss: ',f='(a,$)'
	ptime,tm_in_regauss

	return, outstruct
end



pro admom::memcache_psfield_full, str, index

	common admom_psfield_block_full, run, rerun, camcol, field, psfield, psp

	if n_elements(run) eq 0 then begin
		doread = 1
	endif else begin
		if str[index].run ne run or str[index].rerun ne rerun $
				or str[index].camcol ne camcol $
				or str[index].field ne field then begin
			doread=1
		endif else begin
			doread=0
		endelse
	endelse

	if doread then begin
		run=str[index].run
		rerun=str[index].rerun
		camcol=str[index].camcol
		field=str[index].field

		sf=obj_new('sdss_files')
		if n_elements(psp) ne 0 then begin
			ptr_free, psp
		endif

		; read psf statistics
		psfield = sf->psfield_read(run,camcol,field,$
			rerun=rerun,verbose=0, extension=6)

		; read the kl decompositions
		psp = sf->psfield_read(run,camcol,field,$
			rerun=rerun,verbose=0)

		obj_destroy, sf

	endif

end



pro admom::test_convolve

	; begin with round object
	ixx = 1.0
	ixy = 0.0
	iyy = 1.0

	imsize=[21,21]
	im = mom2gauss(ixx,ixy,iyy,imsize)

	;self->admom, $
	admom, $
		im, $
		10., 10., 0.0, 1.0, 1.0, $
		tixx, tiyy, tixy, momerr, rho4, whyflag

	print,'testing recover of original'
	print,ixx,tixx
	print,ixy,tixy
	print,iyy,tiyy


	kixx = 3.0
	kixy = 1.0
	kiyy = 1.0

	kernel = mom2gauss(kixx, kixy, kiyy, imsize)
	;!p.multi=[0,0,2]
	;rdis,kernel
	;rdis,rotate(kernel,2)
	;!p.multi=0

	;self->admom, $
	admom, $
		kernel, $
		10., 10., 0.0, 1.0, 2.0, $
		tixx, tiyy, tixy, momerr, rho4, whyflag

	print,'testing recover of kernel'
	print,kixx, tixx
	print,kixy, tixy
	print,kiyy, tiyy



	; expected values after convolution
	expixx = ixx + kixx
	expixy = ixy + kixy
	expiyy = iyy + kiyy

	; now test using convol() builtin from IDL
	imconvol = convol(im, kernel, /edge_truncate)

	;self->admom, $
	admom, $
		imconvol, $
		10., 10., 0.0, 1.0, 2.0, $
		tixx, tiyy, tixy, momerr, rho4, whyflag

	print,'testing IDL builtin convol with /edge_truncate'
	print,expixx, tixx
	print,expixy, tixy
	print,expiyy, tiyy



	; now test using convolve(/corr) from goddard
	imconvolve = convolve(im, kernel, /corr)

	;self->admom, $
	admom, $
		imconvolve, $
		10., 10., 0.0, 1.0, 2.0, $
		tixx, tiyy, tixy, momerr, rho4, whyflag

	print,'testing goddard convolve with /corr'
	print,expixx, tixx
	print,expixy, tixy
	print,expiyy, tiyy
	print,'image residual total(imconvolve-imconvol)'
	print,total(imconvolve-imconvol)


	; now test using convolve() from goddard *WITHOUT /CORR*
	imconvolve_nocorr = convolve(im, kernel)

	;self->admom, $
	admom, $
		imconvolve_nocorr, $
		10., 10., 0.0, 1.0, 2.0, $
		tixx, tiyy, tixy, momerr, rho4, whyflag

	print,'testing goddard convolve without /corr'
	print,expixx, tixx
	print,expixy, tixy
	print,expiyy, tiyy
	print,'image residual total(imconvolve_nocorr-imconvol)'
	print,total(imconvolve_nocorr-imconvol)

	; now test using convolve( rotate(kernel,2))
	imconvolve_rot = convolve(im, rotate(kernel,2))

	;self->admom, $
	admom, $
		imconvolve_rot, $
		10., 10., 0.0, 1.0, 2.0, $
		tixx, tiyy, tixy, momerr, rho4, whyflag

	print,'testing goddard convolve with rotate(kernel,2)'
	print,expixx, tixx
	print,expixy, tixy
	print,expiyy, tiyy
	print,'image residual total(imconvolve_rot-imconvol)'
	print,total(imconvolve_rot-imconvol)



end





function admom::make_f0, $
		ixx, ixy, iyy, ixx_psf, ixy_psf, iyy_psf, $
		imsize, cen, counts, $
		status=status

	status=1

	; we want to go to n sigma on each side, so the width will be
	; n*sigma*2 = n*ixx*2 or n*iyy*2 whichever is greater

	ixxf0 = ixx-ixx_psf
	ixyf0 = ixy-ixy_psf
	iyyf0 = iyy-iyy_psf

	det = ixxf0*iyyf0 - ixyf0^2

	; we will cut much harder on this in other places
	if det gt self.detf0_tol then begin
		gauss=mom2gauss(ixxf0, ixyf0, iyyf0, imsize, cen=cen,counts=counts)
		status=0
		return, gauss
	endif else begin
		return, -1
	endelse

end




function admom::make_epsilon, psf, ixx, ixy, iyy, max_size, $
		gauss=gauss
	; make a model image for the fit gaussian and subtract it from the
	; psf.  Make sure the size is within max_size limits.  

	; This becomes a convolution kernel on our simplified model for the galaxy

	; currently this is using a too-large psf in most cases.  We could easily
	; reduce the size to say a 3.5 sigma or 4 sigma region.

	sz=size(psf)
	sx=sz[1]
	sy=sz[2]

	if sx gt max_size[0] then begin
		sx = max_size[0]
	endif
	if sy gt max_size[1] then begin
		sy = max_size[1] 
	endif

	gauss=mom2gauss($
		ixx, ixy, iyy, [sx,sy], counts=1)

	psf_cen = [(sz[1]-1)/2., (sz[2]-1)/2.]
	xrange = [psf_cen[0] - (sx-1.)/2., psf_cen[0]+(sx-1.)/2.]
	yrange = [psf_cen[1] - (sy-1.)/2., psf_cen[1]+(sy-1.)/2.]

	tpsf = psf[xrange[0]:xrange[1], yrange[0]:yrange[1]]

	;new_sz = size(tpsf)
	;if new_sz[0] ne sx or new_sz[1] ne sy then message,'new psf wrong size'

	tpsf = tpsf/total(tpsf)
	epsilon = tpsf - gauss
	return, epsilon

end



pro admom::calculate_corr, str, regauss=regauss, index=index, band=band

	if n_elements(index) eq 0 then index=0
	if n_elements(band) eq 0 then band=0


	if keyword_set(regauss) then begin
		ixx=str[index].ixx_rg[band]
		iyy=str[index].iyy_rg[band]
		ixy=str[index].ixy_rg[band]
		s2o=str[index].s2_rg[band]
		rho4 = (str[index].a4_rg[band]+1)*2
	endif else begin
		ixx=str[index].ixx[band]
		iyy=str[index].iyy[band]
		ixy=str[index].ixy[band]
		s2o=str[index].s2[band]
		rho4 = (str[index].a4[band]+1)*2
	endelse
	ixx_psf=str[index].ixx_psf[band]
	ixy_psf=str[index].ixy_psf[band]
	iyy_psf=str[index].iyy_psf[band]
	s2psf=str[index].s2_psf[band]
	rho4_psf = (str[index].a4_psf[band]+1)*2

	T=ixx+iyy
	e1 = (ixx-iyy)/T
	e2 = 2.0*ixy/T

	Tpsf=ixx_psf+iyy_psf
	e1_psf = (ixx_psf-iyy_psf)/Tpsf
	e2_psf = 2.*ixy_psf/Tpsf

    if 1 then begin
        compea4, $
            e1, e2, T, rho4, $
            e1_psf, e2_psf, Tpsf, rho4_psf, $
            e1out, e2out, Rout, flags

        if keyword_set(regauss) then begin
            str[index].e1_rg[band] = e1out
            str[index].e2_rg[band] = e2out
            str[index].r_rg[band] = Rout
            str[index].whyflag_rg[band] += flags
        endif else begin
            str[index].e1[band] = e1out
            str[index].e2[band] = e2out
            str[index].r[band] = Rout
            str[index].whyflag[band] += flags
        endelse

    endif else begin

	
        R = 1-s2psf/s2o

        if R gt 0 then begin

            e1 = (e1 - (1-R)*e1_psf)/R
            e2 = (e2 - (1-R)*e2_psf)/R

            if keyword_set(regauss) then begin
                str[index].e1_rg[band] = e1
                str[index].e2_rg[band] = e2
                str[index].r_rg[band] = r
            endif else begin
                str[index].e1[band] = e1
                str[index].e2[band] = e2
                str[index].r[band] = r
            endelse
        endif

    endelse
end

pro admom::mom2e1e2, ixx, ixy, iyy, e1, e2
	e1=(ixx-iyy)/(ixx+iyy)	
	e2=2.0*ixy/(ixx+iyy)
end

pro admom::ellip2mom, e, theta, T, ixx, ixy, iyy
	e1 = e*cos(2*theta)
	e2 = e*sin(2*theta)

	ixy = e2*T/2.0

	ixx = (1+e1)*T/2.0
	iyy = (1-e1)*T/2.0
end




function admom::do_regauss, $
		image, cen, sky, skysig, wguess, $
		image_psf, cen_psf, $
		instruct=instruct, index=index, band=band

	common _tmblock, tmconv, tm_in_regauss


	if size(image,/tname) ne 'FLOAT' then message,'Input images must be float'
	if size(image_psf,/tname) ne 'FLOAT' then message,'Input psf must be float'

	if n_elements(instruct) eq 0 then begin
		str = self->regauss_struct(1,/oneband)
		index=0
		band=0
	endif else begin
		; we will start with the input struct and return it, preserving
		; existing information
		if n_elements(index) eq 0 then index=0
		if n_elements(band) eq 0 then band=0

		str = instruct[index]
	endelse

	tm00=systime(1)

	; now the object
	;self->admom, $
	admom, $
		image, cen[0], cen[0], sky, skysig, wguess, $
		ixx, iyy, ixy, momerr, rho4, whyflag

	str.whyflag[band] = whyflag
	if whyflag eq 0 then begin
		a4 = rho4/2-1
		if ixx eq -9999.0 then message,'What?'
		str.ixx[band] = ixx
		str.ixy[band] = ixy
		str.iyy[band] = iyy
		str.momerr[band] = momerr
		str.a4[band] = a4
		str.s2[band] = (ixx+iyy)*(1-a4)/(1+a4)
	endif


	; first do the psf
	sky_psf = 0.0
	skysig_psf = 5.5 ; doesn't matter
	wguess_psf = 1.5
	;self->admom, $
	admom, $
		image_psf, cen_psf[0], cen_psf[0], sky_psf, skysig_psf, wguess_psf, $
		ixx, iyy, ixy, momerr, rho4, whyflag

	str.whyflag_psf[band] = whyflag
	if whyflag eq 0 then begin
		a4 = rho4/2-1
		str.ixx_psf[band] = ixx
		str.ixy_psf[band] = ixy
		str.iyy_psf[band] = iyy
		str.a4_psf[band] = a4
		str.s2_psf[band] = (ixx+iyy)*(1-a4)/(1+a4)
	endif



	; to proceed we require good measuremetns for both the PSF and the object
	if str.whyflag[band] eq 0 and str.whyflag_psf[band] eq 0 then begin

		; apply the ordinary correction just for comparisons later
        ; these are stored in e1, e2, etc. instead of e1_rg, etc.
		self->calculate_corr, str, band=band

		; construct the f0 parameters
		str.ixx_f0[band] = str.ixx[band] - str.ixx_psf[band]
		str.ixy_f0[band] = str.ixy[band] - str.ixy_psf[band]
		str.iyy_f0[band] = str.iyy[band] - str.iyy_psf[band]

		str.detf0[band] = $
			str.ixx_f0[band]*str.iyy_f0[band] - str.ixy_f0[band]^2

		; make sure the new covariance matrix is positive definite
		; make sure the sizes in each direction are positive
		; make sure the resolution parameter is sane
		if str.detf0[band] gt self.detf0_lim $
				and str.ixx_f0[band] gt 0.01 $
				and str.iyy_f0[band] gt 0.01 $
				and str.r[band] gt 0.1  then begin


			; now develop the re-gaussianized image.  f0 is created on the
			; same image area as the atlas image
			imsize=size(image)
			imsize=imsize[1:2]

			imcounts=total(image)
			f0 = self->make_f0($
				str.ixx[band], str.ixy[band], str.iyy[band], $
				str.ixx_psf[band], str.ixy_psf[band], str.iyy_psf[band], $
				imsize,cen,imcounts)
			; Subtract our best gaussian fit from the adaptive moments from
			; the psf
			epsilon = self->make_epsilon($
				image_psf,$
				str.ixx_psf[band], str.ixy_psf[band], str.iyy_psf[band], $
				imsize)

			tm1=systime(1)
			f0conv = convolve(f0, epsilon)
			;f0conv = convol(f0, epsilon)
			tmconv += systime(1)-tm1

			iprime = image-f0conv

			wguess_regauss = str.ixx[band]

			;self->admom, $
			admom, $
				iprime, cen[0], cen[1], sky, skysig, wguess_regauss, $
				ixx, iyy, ixy, momerr, rho4, whyflag

			str.whyflag_rg[band] = whyflag

			if whyflag eq 0 then begin
				a4 = rho4/2-1
				str.ixx_rg[band] = ixx
				str.ixy_rg[band] = ixy
				str.iyy_rg[band] = iyy
				str.momerr_rg[band] = momerr
				str.a4_rg[band] = a4
				str.s2_rg[band] = (ixx+iyy)*(1-a4)/(1+a4)

				; regauss correction.  This will set e*_rg, r_rg and possibly
                ; additional flags in whyflag_rg
				self->calculate_corr, str, band=band, /regauss

				if 0 then begin
					print, $
						index, band, $
						str.e1[band],str.e2[band],str.e1_rg[band],str.e2_rg[band]
					stop
				endif
			endif

		endif
	endif

	tm_in_regauss+=systime(1)-tm00
	return, str
end
pro admom::do_regauss_atlas, str, outstruct, index, band, astr, psf


	common admom_psfield_block_full, run, rerun, camcol, field, psfield, psp

	;; get inputs
	case band of
		0: atlas = float(astr.imu)
		1: atlas = float(astr.img)
		2: atlas = float(astr.imr)
		3: atlas = float(astr.imi)
		4: atlas = float(astr.imz)
		else: message,'bad band: '+strn(band)
	endcase

	xcen = str[index].colc[band] - astr.col0[band] - 0.5
	ycen = str[index].rowc[band] - astr.row0[band] - 0.5
	cen = [xcen,ycen]
	sky = 0.0
	skysig = psfield.skysig[band]

	if tag_exist(str,'m_rr_cc') then begin
		; bound the guesses, even thogh there are real objects that
		; might have ixx+iyy greater than 20
		wguess = (str[index].m_rr_cc[band]/2.0) > 1.0 < 20.0
	endif else begin
		wguess = 2.0
	endelse

	psfsize=(size(psf))[1:2]
	cen_psf = [(psfsize[0]-1.0)/2.0, (psfsize[1]-1.0)/2.0]

	outstruct[index] = self->do_regauss( $
		atlas, cen, sky, skysig, wguess, $
		psf, cen_psf, $
		instruct=outstruct[index], band=band)

	return

end



function admom::regauss_struct, n, oneband=oneband

	if keyword_set(oneband) then begin
		defval=-9999.0
		defval2 = 9999.0
		whyval=2L^0
	endif else begin
		defval=replicate(-9999.0,5)
		defval2 = replicate(9999.0, 5)
		whyval=replicate(2L^0,5)
	endelse
	outstruct = { $
		run: 0, $
		rerun: '', $
		camcol: 0, $
		field: 0, $
		id: 0L, $
		ixx: defval, $
		iyy: defval, $
		ixy: defval, $
		a4: defval, $
		s2: defval, $
		momerr: defval2, $
		whyflag: whyval, $
		$
		ixx_psf: defval, $
		iyy_psf: defval, $
		ixy_psf: defval, $
		a4_psf: defval, $
		s2_psf: defval, $
		whyflag_psf: whyval, $
		$
		ixx_rg: defval, $
		iyy_rg: defval, $
		ixy_rg: defval, $
		a4_rg: defval, $
		s2_rg: defval, $
		momerr_rg: defval2, $
		whyflag_rg: whyval, $
		$
		ixx_f0: defval, $
		ixy_f0: defval, $
		iyy_f0: defval, $
		detf0: defval, $
		$
		r: defval, $
		e1: defval, $
		e2: defval, $
		$
		r_rg: defval, $
		e1_rg: defval, $
		e2_rg: defval $
	}
	return, replicate(outstruct, n)
end






pro admom::get_test_data, str, rg
	if n_elements(str) eq 0 then begin
		columns=['run','rerun','camcol','field','id',$
			'colc','rowc',$
			'm_rr_cc','objc_type','modelflux','extinction']
		f='~/tmp/all-000756-2-301.fits'
		str=mrdfits(f,1,columns=columns)
	endif
	if n_elements(rg) eq 0 then begin
		rg=mrdfits('~/tmp/regauss-000756-2-301.fits',1)
	endif
end






pro admom::explore_regauss, str=str, rg=rg, wgalbig=wgalbig

	self->get_test_data, str, rg

	rmag=22.5-2.5*alog10(str.modelflux[2] > 0.001) - str.extinction[2]
	rmaglim=22.0
	w=where(rmag lt rmaglim and rg.whyflag[2] eq 0 and rg.whyflag_psf[2] eq 0, nw)
	wgal = where(str[w].objc_type eq 3)
	wgal = w[wgal]
	wgalbig = where(rg[wgal].r[2] gt 1.0/3.0, nbig)
	wgalbig = wgal[wgalbig]

	bdet1=0.0
	wgalbig_baddet = where(rg[wgalbig].detf0[2] lt bdet1, nbad_det)
	bdet2=0.1
	wgalbig_baddet2 = where(rg[wgalbig].detf0[2] lt bdet2, nbad_det2)
	bdet3=0.5
	wgalbig_baddet3 = where(rg[wgalbig].detf0[2] lt bdet3, nbad_det3)
	bdet4=1.0
	wgalbig_baddet4 = where(rg[wgalbig].detf0[2] lt bdet4, nbad_det4)
	print,'fraction of big with bad det < ',bdet1,float(nbad_det)/nbig
	print,'fraction of big with bad det < ',bdet2,float(nbad_det2)/nbig
	print,'fraction of big with bad det < ',bdet3,float(nbad_det3)/nbig
	print,'fraction of big with bad det < ',bdet4,float(nbad_det4)/nbig

	begplot,'~/tmp/test-regauss.ps',/color

	!p.multi=[0,0,2]
	rmin=0.0
	rmax=1.0
	plothist, rg[w].r[2],min=rmin,max=rmax,bin=0.01, xtitle='R'
	plothist, rg[wgal].r[2], min=rmin, max=rmax, bin=0.01, $
		/over, color='red'
	plothist, rg[wgalbig].r[2], min=rmin, max=rmax, bin=0.01, $
		/over, color='blue'

	dmin=-2.0
	dmax=25.0
	plothist, rg[w].detf0[2], min=dmin, max=dmax, bin=0.1,/ylog, $
		yrange=[80,2.e4], ystyle=3, xtitle='det(f0)'
	plothist, rg[wgal].detf0[2], min=dmin, max=dmax, bin=0.1, $
		/over, color='red';, /norm
	plothist, rg[wgalbig].detf0[2], min=dmin, max=dmax, bin=0.1, $
		/over, color='blue';, /norm

	!p.multi=0


	endplot



	detlim=0.2
	wplot=where(rg[wgalbig].detf0[2] gt detlim and rg[wgalbig].detf0[2] lt detlim+0.1,nplot)
	wplot=wgalbig[wplot]

	ixx=rg[wplot].ixx[2]-rg[wplot].ixx_psf[2]
	ixy=rg[wplot].ixy[2]-rg[wplot].ixy_psf[2]
	iyy=rg[wplot].iyy[2]-rg[wplot].iyy_psf[2]

	begplot,'~/tmp/test-detf0-'+string(detlim,f='(f0.1)')+'.ps',/color
	!p.multi=[0,0,2]
	nplot = nplot < 100
	for i=0L, nplot-1 do begin
		g=mom2gauss(ixx[i],ixy[i],iyy[i],[31,31])
		loadct,0
		rdis,g,/scale

		simpctable
		pplot,g[15,*],psym=8
		pplot,/over,g[15,*],color='red'
		pplot,/over,g[*,15],psym=8
		pplot,/over,g[*,15],color='blue'
	endfor
	endplot
	!p.multi=0

end









pro admom__define

  struct = {$
             admom, $
			 detf0_lim: 0.0, $
			 detf0_tol: 0.0 $
           }

end 
