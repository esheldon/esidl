function regauss_sim::init
    return,self->regauss::init()
end


;
;
;  Simulations
;
;

; note: for exp models, the ellip in the name may not match
; the ellip in the etrue field but it's in the ballpark
pro regauss_sim::run, psfmodel, galmodel, sizerat2

    self->model_check, psfmodel, galmodel

    dir=self->dir(psfmodel)
    if not file_test(dir) then file_mkdir, dir

	; sizerat2 is the ratio of (psfsize/galsize)^2

	; ellipticities to try
	;ellipvals= [0.1,0.2,0.3,0.4,0.5,0.6,0.7]
	ellipvals=arrscl(findgen(20), 0.05, 0.8)

	nell = n_elements(ellipvals)

	; number of times we will realize each ellipticity
	; each one gets a different psf
	ntrial = 100LL

	; number of random orientations that an object will get
    ; for each trial psf
	nrand=100LL

	for i=0LL, nell-1 do begin
		ellip = ellipvals[i]
		print,'ellip: ',ellip

        outfile=self->outfile(psfmodel, galmodel, sizerat2, ellip)
		print,'will write to file: ',outfile

        j=0LL
        while (j lt ntrial) do begin
            if psfmodel eq 'sdss' then begin
                tres = self->sdss_psf_realizations(galmodel, sizerat2, ellip, nrand=nrand, $
                    status=status)
            endif else begin
                ;tres = self->gauss_psf_realizations(galmodel, sizerat2, ellip, nrand=nrand, $
                ;    status=status)
                tres = self->sim_psf_realizations(psfmodel, galmodel, sizerat2, ellip, nrand=nrand, $
                    status=status)
            endelse

            if status eq 0 then begin
                if j eq 0 then append=0
                write_idlstruct, tres, outfile, append=append, status=status
                if status ne 0 then message,'Failed to write file'
                append=1
                print,'.',f='($,a)'
                j += 1
            endif

        endwhile
        print
	endfor
    print,'Done'

end

function regauss_sim::dir, psfmodel
	dir=getenv('REGAUSSIM_DIR')
    dir = filepath(root=dir, psfmodel+'-psf')
    return, dir
end
function regauss_sim::outfile, psfmodel, galmodel, sizerat2, ellip
    dir = self->dir(psfmodel)
    if size(ellip, /tname) eq 'STRING' then begin
        ellipstr = ellip
    endif else begin
        ellipstr = string(ellip, f='(f0.2)')
    endelse
    file=string($
        f='("regauss-sim-",a,"-",a,"-sizerat2-",f0.2,"-ell-",a,".st")',$
             psfmodel,galmodel,sizerat2,ellipstr)
    file=path_join(dir, file)
    return,file
end

pro regauss_sim::model_check, psfmodel, galmodel
    if not in(['gauss','exp'], galmodel) then begin
        message,'galmodel must be gauss or exp'
    endif
    if not in(['gauss','dgauss','sdss'], psfmodel) then begin
        message,'psf model must be gauss or sdss'
    endif
end

function regauss_sim::sdss_psf_realizations, galmodel, sizerat2, ellip, index=index, nrand=nrand, $
        status=status

    common _regauss_sim_sdss_psf_random_block, str
	common admom_psfield_block_full, run, rerun, camcol, field, psfield, psp

    if n_elements(str) eq 0 then begin
        str = self->random_objects_read()
    endif
    ; some times when the sizerat2 is big we fail and need to 
    ; draw another psf.  Status will be 1 then
    status=1

    self->model_check, 'sdss', galmodel

    nrand=long64(nrand)

    ; for a random object, get it's PSF info and then
    ; choose nrand ellipticity orientations

	; sizerat2 is the ratio of (psfsize/galsize)^2
	; str is a photo structure so we can get a realistic PSF
	; ellip is the total ellipticity of the galaxy to create


	if n_elements(nrand) eq 0 then nrand=1LL



	; get a realistic PSF for a random object
	nstr=n_elements(str)
	if n_elements(index) eq 0 then begin
		; pick a random object to get the psf info for
		index = long( randomu(seed)*(nstr-1) )
	endif

	;print,'using index=',index
	self->regauss::memcache_psfield_full, str, index, verbose=0
	band = 2
	psf = sdss_psfrec($
		*psp[band], $
		str[index].rowc[band], str[index].colc[band],$
		counts=1.0, /trim)
	psfsize=(size(psf))[1:2]
	cenpsf = [(psfsize[0]-1)/2., (psfsize[1]-1)/2.]

	; measure the psf
	admom, $
		psf, cenpsf[0], cenpsf[1], 0.0, 5.5, 1.5, $
		ixx_psf, iyy_psf, ixy_psf, momerr_psf, rho4_psf, whyflag_psf

	Tpsf = (ixx_psf + iyy_psf)

	Tgal = Tpsf/sizerat2


	; now make the random orientations


	struct = {$
		sizerat2: sizerat2, $
        sizerat2_actual: -9999.0, $
		e1true: -9999.0, $
		e2true: -9999.0, $
		etrue: -9999.0, $
		e1meas: -9999.0, $
		e2meas: -9999.0, $
		emeas: -9999.0, $
        whyflag: -9999L $
	}
	struct = replicate(struct, nrand)


	imsize=[41,41]
	cen=[ (imsize[0]-1)/2., (imsize[1]-1)/2. ]
	sky=0.0
	skysig=5.5 ; does not matter

    ; loop until we get enough successes
    i=0LL
    ntry=0LL
    trymax = nrand*100
    while (i lt nrand) and (ntry lt trymax) do begin
    ;while (i lt nrand) do begin
	;for i=0L, nrand-1 do begin

        ntry += 1

	    theta = randomu(seed)*2.0*!pi
        self->regauss::ellip2mom, ellip, theta, Tgal, ixx, ixy, iyy
        object = mom2disk(galmodel,ixx, ixy, iyy, imsize, counts=1.0)

        if 0 then begin
            case galmodel of
                'gauss': begin
                    e1true = (ixx-iyy)/(ixx+iyy)
                    e2true = 2.*ixy/(ixx+iyy)
                    etrue = ellip
                    ; so we pass the check below
                    twhyflag=0
                end
                'exp': begin

                    wguess=Tgal/2.0
                    admom, $
                        object, cen[0], cen[1], sky, skysig, wguess, $
                        tixx, tiyy, tixy, tmomerr, trho4, twhyflag
                    if twhyflag ne 0 then begin
                        ;print,'admom failed on exp object'
                    endif
                    tT = tixx+tiyy
                    e1true = (tixx-tiyy)/tT
                    e2true = 2.0*tixy/tT
                    etrue = sqrt( e1true^2 + e2true^2 )
                end
            endcase
        endif

        ; use truth as what we get from adaptive moments
        ; before convolution.  This avoids any other
        ; artifacts from pixelization, etc. confusing
        ; our attempt to correct for the PSF
        wguess=Tgal/2.0
        admom, $
            object, cen[0], cen[1], sky, skysig, wguess, $
            tixx, tiyy, tixy, tmomerr, trho4, twhyflag


        ; only use converged values; otherwise take another
        ; trial
        if twhyflag eq 0 then begin

            tT = tixx+tiyy
            e1true = (tixx-tiyy)/tT
            e2true = 2.0*tixy/tT
            etrue = sqrt( e1true^2 + e2true^2 )


            ; now convolve the object with the psf
            object_convolved = convolve(object, psf)

            ;wguess = ixx
            wguess = Tgal/2.0
            res = self->regauss::do_regauss($
                object_convolved, cen, sky, skysig, wguess, $
                psf, cenpsf)

            if (res.whyflag eq 0) then begin
                struct[i].sizerat2_actual = Tpsf/tT
                struct[i].e1true = e1true
                struct[i].e2true = e2true
                struct[i].etrue = etrue

                struct[i].e1meas = res.e1_rg
                struct[i].e2meas = res.e2_rg
                struct[i].emeas = sqrt(res.e1_rg^2 + res.e2_rg^2)
                struct[i].whyflag = res.whyflag
                i += 1
            endif
        endif
    endwhile

    if ntry ge trymax then begin
        ;print,'  reached max tries',trymax
    endif else begin
        status=0
    endelse

    return, struct
end

function regauss_sim::sim_psf_realizations, psfmodel, galmodel, sizerat2, ellip, nrand=nrand, status=status, $
        psf_fwhm=psf_fwhm, $
        doplot=doplot

    ; sizerat2 is (psfsize/galsize)^2
    ; nrand is the number of random angles to realize

    ; some times when the sizerat2 is big we fail and need to 
    ; draw another psf.  Status will be 1 then
    status=1
    self->model_check, psfmodel, galmodel

	if n_elements(nrand) eq 0 then nrand=1LL else nrand=long64(nrand)

    ; take the psf size and ellip to be fixed for now
    fac = 2*sqrt(2*alog(2))
    ; arcsec per pixel
    pixscale = 0.4

    ; arcsec
    ;psf_fwhm = 1.3
    ; remove pixelization issues
    if n_elements(psf_fwhm) eq 0 then psf_fwhm = 2.5

    ; sigma in pixels
    psf_sigma = psf_fwhm/fac/pixscale

    psf_ellip = 0.1


    ;sigma = FWHM/fac
    ;mom = 2*(fwhm/fac)^2 = 2*sigma^2
    ; thus sigma = 

    ; T is ixx+iyy
    Tpsf_input = fwhm2mom(psf_fwhm, pixscale=pixscale)

	psf_theta = randomu(seed)*2.0*!pi

    if psfmodel eq 'gauss' then begin
        self->regauss::ellip2mom, psf_ellip, psf_theta, Tpsf_input, psf_ixx_input, psf_ixy_input, psf_iyy_input

        psf_imsize = [2*4*psf_sigma, 2*4*psf_sigma]
        psf = mom2disk('gauss', $
                       psf_ixx_input, psf_ixy_input, psf_iyy_input, $
                       psf_imsize, counts=1.0, cen=psf_cen)
    endif else begin
        ; twice the sigma, 4 times the T
        Tpsf_input2 = 4*Tpsf_input
        ; 80% of flux in the first
        psf = self->regauss::double_gauss(psf_ellip, psf_theta, Tpsf_input, Tpsf_input2, 0.8)
        psf_imsize = size(psf, /dim)
        psf_cen = [0.0, 0.0] + (psf_imsize[0]-1.0)/2.0
    endelse

    ; now get the actual measurements.  Since this is a gaussian, it should match extremely
    ; well

	; measure the psf
	admom, $
		psf, psf_cen[0], psf_cen[1], 0.0, 5.5, 1.5, $
		psf_ixx, psf_iyy, psf_ixy, psf_momerr, psf_rho4, psf_whyflag

	Tpsf = (psf_ixx + psf_iyy)
	Tgal = Tpsf/sizerat2

    gal_fwhm = mom2fwhm(Tgal, pixscale=pixscale) ; in arcsec
    gal_sigma = gal_fwhm/fac/pixscale ; back in pixels


	; now make the random galaxy orientations

	struct = {$
		sizerat2: sizerat2, $
        sizerat2_actual: -9999.0, $
		e1true: -9999.0, $
		e2true: -9999.0, $
		etrue: -9999.0, $
        e1convolved: -9999.0, $
        e2convolved: -9999.0, $
        econvolved: -9999.0, $
		e1meas: -9999.0, $
		e2meas: -9999.0, $
		emeas: -9999.0, $
        whyflag: -9999L, $
		e1meas_rg: -9999.0, $
		e2meas_rg: -9999.0, $
		emeas_rg: -9999.0, $
        whyflag_rg: -9999L $
	}
	struct = replicate(struct, nrand)

    ; note psf_size above is [11,11] for fwhm=1.3
	;imsize=[41,41]
    if galmodel eq 'exp' then begin
        sigfac = 5.0
    endif else begin
        sigfac = 4.0
    endelse
    imsize = 2*sigfac*[gal_sigma, gal_sigma]
    ;if imsize[0] lt psf_imsize[0] then begin
    ;    imsize = psf_imsize
    ;endif
	cen=[ (imsize[0]-1)/2., (imsize[1]-1)/2. ]
	sky=0.0
	skysig=5.5 ; does not matter

    ;if psf_imsize[0] gt imsize[0] then message,'psf is larger than galaxy'

    ; loop until we get enough successes
    i=0LL
    ntry=0LL
    trymax = nrand*100
    while (i lt nrand) and (ntry lt trymax) do begin
    ;while (i lt nrand) do begin
	;for i=0L, nrand-1 do begin

        ntry += 1

	    theta = randomu(seed)*2.0*!pi
        self->regauss::ellip2mom, ellip, theta, Tgal, ixx, ixy, iyy
        object = mom2disk(galmodel,ixx, ixy, iyy, imsize, counts=1.0)


        if 0 then begin
            case galmodel of
                'gauss': begin
                    e1true = (ixx-iyy)/(ixx+iyy)
                    e2true = 2.*ixy/(ixx+iyy)
                    etrue = ellip
                    ; so we pass the check below
                    twhyflag=0
                end
                'exp': begin

                    wguess=Tgal/2.0
                    admom, $
                        object, cen[0], cen[1], sky, skysig, wguess, $
                        tixx, tiyy, tixy, tmomerr, trho4, twhyflag
                    if twhyflag ne 0 then begin
                        ;print,'admom failed on exp object'
                    endif
                    tT = tixx+tiyy
                    e1true = (tixx-tiyy)/tT
                    e2true = 2.0*tixy/tT
                    etrue = sqrt( e1true^2 + e2true^2 )
                end
            endcase
        endif

        ; use truth as what we get from adaptive moments
        ; before convolution.  This avoids any other
        ; artifacts from pixelization, etc. confusing
        ; our attempt to correct for the PSF
        wguess=Tgal/2.0
        admom, $
            object, cen[0], cen[1], sky, skysig, wguess, $
            tixx, tiyy, tixy, tmomerr, trho4, twhyflag


        ; only use converged values; otherwise take another
        ; trial
        if twhyflag eq 0 then begin

            tT = tixx+tiyy
            e1true = (tixx-tiyy)/tT
            e2true = 2.0*tixy/tT
            etrue = sqrt( e1true^2 + e2true^2 )


            ; now convolve the object with the psf
            object_convolved = convolve(object, psf)

            ;wguess = ixx
            wguess = Tgal/2.0
            ;res = self->regauss::do_regauss($
            ;    object_convolved, cen, sky, skysig, wguess, $
            ;    psf, cenpsf)
            res = self->regauss::process_image($
                object_convolved, cen, sky, skysig, wguess, $
                psf, psf_cen)

            if (res.whyflag eq 0) then begin
                struct[i].sizerat2_actual = Tpsf/tT
                struct[i].e1true = e1true
                struct[i].e2true = e2true
                struct[i].etrue = etrue

                struct[i].e1convolved = (res.ixx - res.iyy)/(res.ixx + res.iyy)
                struct[i].e2convolved = 2*res.ixy/(res.ixx + res.iyy)
                struct[i].econvolved = sqrt(struct[i].e1convolved^2 + struct[i].e2convolved^2)

                struct[i].e1meas = res.e1
                struct[i].e2meas = res.e2
                struct[i].emeas = sqrt(res.e1^2 + res.e2^2)
                struct[i].whyflag = res.whyflag

                struct[i].e1meas_rg = res.e1_rg
                struct[i].e2meas_rg = res.e2_rg
                struct[i].emeas_rg = sqrt(res.e1_rg^2 + res.e2_rg^2)
                struct[i].whyflag_rg = res.whyflag_rg
                i += 1
            endif
        endif
    endwhile

    if ntry ge trymax then begin
        print,'  reached max tries',trymax
    endif else begin
        status=0

        if keyword_set(doplot) then begin
            !p.multi=[0,2,0]
            !p.charsize=2

            w=where(struct.whyflag eq 0, nw)
            print,'Found ',nw,' with good whyflag',format='(a,i0,a)'
            if nw gt 2 then begin
                e1diff = struct[w].e1meas - struct[w].e1true
                e2diff = struct[w].e2meas - struct[w].e2true
                std1 = stdev(e1diff)
                std2 = stdev(e2diff)
                std = max([std1,std2])
                binsize = std/8.0

                plothist, e1diff, binsize=binsize, ystyle=3
                plothist, e2diff, /over, binsize=binsize, color='green'
                plegend,'CH'
                plegend, ['e1-e1true','e2-e2true'], line=0, color=[!p.color, c2i('green')], $
                    /right
            endif

            w=where(struct.whyflag_rg eq 0, nw)
            print,'Found ',nw,' with good whyflag_rg',format='(a,i0,a)'
            if nw gt 2 then begin
                e1diff = struct[w].e1meas_rg - struct[w].e1true
                e2diff = struct[w].e2meas_rg - struct[w].e2true
                std1 = stdev(e1diff)
                std2 = stdev(e2diff)
                std = max([std1,std2])
                binsize = std/8.0

                plothist, e1diff, binsize=binsize, ystyle=3
                plothist, e2diff, /over, binsize=binsize, color='green'

                plegend,'regauss'
                plegend, ['e1-e1true','e2-e2true'], line=0, color=[!p.color, c2i('green')], $
                    /right
            endif



            !p.multi=0
            !p.charsize=1
        endif

    endelse

    return, struct
end



function regauss_sim::gauss_psf_realizations, galmodel, sizerat2, ellip, nrand=nrand, status=status, $
        doplot=doplot

    psfmodel = 'gauss'
    ; some times when the sizerat2 is big we fail and need to 
    ; draw another psf.  Status will be 1 then
    status=1
    self->model_check, psfmodel, galmodel

	if n_elements(nrand) eq 0 then nrand=1LL else nrand=long64(nrand)

    ; take the psf size and ellip to be fixed for now
    fac = 2*sqrt(2*alog(2))
    ; arcsec per pixel
    pixscale = 0.4

    ; arcsec
    ;psf_fwhm = 1.3
    ; remove pixelization issues
    psf_fwhm = 2.5

    ; sigma in pixels
    psf_sigma = psf_fwhm/fac/pixscale

    psf_ellip = 0.1


    ;sigma = FWHM/fac
    ;mom = 2*(fwhm/fac)^2 = 2*sigma^2
    ; thus sigma = 

    ; T is ixx+iyy
    Tpsf_input = fwhm2mom(psf_fwhm, pixscale=pixscale)

	psf_theta = randomu(seed)*2.0*!pi
    self->regauss::ellip2mom, psf_ellip, psf_theta, Tpsf_input, psf_ixx_input, psf_ixy_input, psf_iyy_input

    psf_imsize = [2*4*psf_sigma, 2*4*psf_sigma]
    psf = mom2disk('gauss', $
                   psf_ixx_input, psf_ixy_input, psf_iyy_input, $
                   psf_imsize, counts=1.0, cen=psf_cen)

    ; now get the actual measurements.  Since this is a gaussian, it should match extremely
    ; well

	; measure the psf
	admom, $
		psf, psf_cen[0], psf_cen[1], 0.0, 5.5, 1.5, $
		psf_ixx, psf_iyy, psf_ixy, psf_momerr, psf_rho4, psf_whyflag

	Tpsf = (psf_ixx + psf_iyy)
	Tgal = Tpsf/sizerat2

    gal_fwhm = mom2fwhm(Tgal, pixscale=pixscale) ; in arcsec
    gal_sigma = gal_fwhm/fac/pixscale ; back in pixels


    if keyword_set(doplot) then begin
        print, 'psf_imsize: ',psf_imsize
        help, Tpsf_input, Tpsf
        help, psf_ixx_input, psf_ixx
        help, psf_ixy_input, psf_ixy
        help, psf_iyy_input, psf_iyy

        !p.multi=[0,2,2]
        tvim2, psf
        legend,psfmodel+' PSF',/right
        plot, psf[*,psf_cen[1]]
        oplot, psf[psf_cen[0],*], color=c2i('green')
        legend,psfmodel+' PSF',/right
    endif

	; now make the random galaxy orientations

	struct = {$
		sizerat2: sizerat2, $
        sizerat2_actual: -9999.0, $
		e1true: -9999.0, $
		e2true: -9999.0, $
		etrue: -9999.0, $
        e1convolved: -9999.0, $
        e2convolved: -9999.0, $
        econvolved: -9999.0, $
		e1meas: -9999.0, $
		e2meas: -9999.0, $
		emeas: -9999.0, $
        whyflag: -9999L, $
		e1meas_rg: -9999.0, $
		e2meas_rg: -9999.0, $
		emeas_rg: -9999.0, $
        whyflag_rg: -9999L $
	}
	struct = replicate(struct, nrand)

    ; note psf_size above is [11,11] for fwhm=1.3
	;imsize=[41,41]
    if galmodel eq 'exp' then begin
        sigfac = 5.0
    endif else begin
        sigfac = 4.0
    endelse
    imsize = 2*sigfac*[gal_sigma, gal_sigma]
    ;if imsize[0] lt psf_imsize[0] then begin
    ;    imsize = psf_imsize
    ;endif
	cen=[ (imsize[0]-1)/2., (imsize[1]-1)/2. ]
	sky=0.0
	skysig=5.5 ; does not matter

    ;if psf_imsize[0] gt imsize[0] then message,'psf is larger than galaxy'

    ; loop until we get enough successes
    i=0LL
    ntry=0LL
    trymax = nrand*100
    while (i lt nrand) and (ntry lt trymax) do begin
    ;while (i lt nrand) do begin
	;for i=0L, nrand-1 do begin

        ntry += 1

	    theta = randomu(seed)*2.0*!pi
        self->regauss::ellip2mom, ellip, theta, Tgal, ixx, ixy, iyy
        object = mom2disk(galmodel,ixx, ixy, iyy, imsize, counts=1.0)

        if i eq 0 and keyword_set(doplot) then begin
            print,'object image size: ',imsize
            tvim2, object
            legend,galmodel+' gal',/right
            plot, object[*,cen[1]]
            oplot, object[cen[0],*], color=c2i('green')
            legend,galmodel+' gal',/right
        endif

        if 0 then begin
            case galmodel of
                'gauss': begin
                    e1true = (ixx-iyy)/(ixx+iyy)
                    e2true = 2.*ixy/(ixx+iyy)
                    etrue = ellip
                    ; so we pass the check below
                    twhyflag=0
                end
                'exp': begin

                    wguess=Tgal/2.0
                    admom, $
                        object, cen[0], cen[1], sky, skysig, wguess, $
                        tixx, tiyy, tixy, tmomerr, trho4, twhyflag
                    if twhyflag ne 0 then begin
                        ;print,'admom failed on exp object'
                    endif
                    tT = tixx+tiyy
                    e1true = (tixx-tiyy)/tT
                    e2true = 2.0*tixy/tT
                    etrue = sqrt( e1true^2 + e2true^2 )
                end
            endcase
        endif

        ; use truth as what we get from adaptive moments
        ; before convolution.  This avoids any other
        ; artifacts from pixelization, etc. confusing
        ; our attempt to correct for the PSF
        wguess=Tgal/2.0
        admom, $
            object, cen[0], cen[1], sky, skysig, wguess, $
            tixx, tiyy, tixy, tmomerr, trho4, twhyflag


        ; only use converged values; otherwise take another
        ; trial
        if twhyflag eq 0 then begin

            tT = tixx+tiyy
            e1true = (tixx-tiyy)/tT
            e2true = 2.0*tixy/tT
            etrue = sqrt( e1true^2 + e2true^2 )


            ; now convolve the object with the psf
            object_convolved = convolve(object, psf)

            ;wguess = ixx
            wguess = Tgal/2.0
            ;res = self->regauss::do_regauss($
            ;    object_convolved, cen, sky, skysig, wguess, $
            ;    psf, cenpsf)
            res = self->regauss::process_image($
                object_convolved, cen, sky, skysig, wguess, $
                psf, psf_cen)

            if (res.whyflag eq 0) then begin
                struct[i].sizerat2_actual = Tpsf/tT
                struct[i].e1true = e1true
                struct[i].e2true = e2true
                struct[i].etrue = etrue

                struct[i].e1convolved = (res.ixx - res.iyy)/(res.ixx + res.iyy)
                struct[i].e2convolved = 2*res.ixy/(res.ixx + res.iyy)
                struct[i].econvolved = sqrt(struct[i].e1convolved^2 + struct[i].e2convolved^2)

                struct[i].e1meas = res.e1
                struct[i].e2meas = res.e2
                struct[i].emeas = sqrt(res.e1^2 + res.e2^2)
                struct[i].whyflag = res.whyflag

                struct[i].e1meas_rg = res.e1_rg
                struct[i].e2meas_rg = res.e2_rg
                struct[i].emeas_rg = sqrt(res.e1_rg^2 + res.e2_rg^2)
                struct[i].whyflag_rg = res.whyflag_rg
                i += 1
            endif
        endif
    endwhile

    if ntry ge trymax then begin
        print,'  reached max tries',trymax
    endif else begin
        status=0
    endelse

    return, struct
end




; read in all the files for this model and approximate
; sizerat2 and return some stats
function regauss_sim::getstats, psfmodel, galmodel, sizerat2

    if not in(['gauss','exp'], galmodel) then begin
        message,'galmodel must be gauss or exp'
    endif
	dir=getenv('REGAUSSIM_DIR')

    pattern = self->outfile(psfmodel, galmodel, sizerat2, '*')

	print,'looking for pattern: ',pattern
	f=file_search(pattern)
	nf=n_elements(f)

	st=replicate({$
        ellip:-9999.0, $
        pdiff:-9999.0, $
        pdiff_sdev:-9999.0,$
        pdiff_rg:-9999.0, $
        pdiff_rg_sdev:-9999.0,$
        sizerat2_actual:-9999.0}, nf)
	for i=0L, nf-1 do begin
        print,f[i]
		t=read_idlstruct(f[i])

        ; whyflag is only -9999 if we hit trymax
        w=where(t.whyflag eq 0, nw)
        if nw gt 0 then begin

            pdiff = (t[w].emeas-t[w].etrue)/t[w].etrue

            sigma_clip, pdiff, mn, sig
            st[i].ellip = t[0].etrue
            st[i].pdiff = mn
            ;st[i].pdiff = median(pdiff)
            st[i].pdiff_sdev = sig

            st[i].sizerat2_actual = median(t[w].sizerat2_actual)

            w2 = where(t[w].whyflag_rg eq 0, n2)
            if n2 gt 0 then begin
                w2=w[w2]
                pdiff_rg = (t[w2].emeas_rg-t[w2].etrue)/t[w2].etrue

                sigma_clip, pdiff_rg, mn, sig
                st[i].pdiff_rg = mn
                ;st[i].pdiff = median(pdiff)
                st[i].pdiff_rg_sdev = sig

            endif
        endif
	endfor

    return, st

end

; plot all the size ratios on a single plot
pro regauss_sim::plot_allone, psfmodel, galmodel, dops=dops
    ;sr2=self->example_sizerat2(/less)
    sr2=self->example_sizerat2()
    nsr2=n_elements(sr2)


    if psfmodel eq 'gauss' then begin
        yrange=[-0.03,0.05]
    endif else if psfmodel eq 'dgauss' then begin
        yrange = [-0.05,0.06]
    endif else begin
        yrange=[-0.01,0.05]
    endelse
    xrange = [0,0.8]

    dir = self->dir(psfmodel)
    psfile=string($
        f='("regauss-sim-",a,"-",a,".eps")',psfmodel,galmodel)
    psfile=filepath(root=dir,psfile)
    if keyword_set(dops) then begin
        begplot,psfile,/encap,/color;, xsize=8.5, ysize=8.5
    endif else begin
        !p.charsize=2
    endelse

    !p.multi=[0,0,2]
    colors=make_rainbow(nsr2)

    plot, [0], /nodata, $
        yrange=yrange, ystyle=3, xrange=xrange, xstyle=3, $
        xtitle='intrinsic ellip', ytitle=textoidl('\deltae/e');, $
        ;aspect=esheldon_config('gratio')

    for i=0L, nsr2-1 do begin

        t=self->getstats(psfmodel, galmodel, sr2[i])
        w=where(t.ellip gt 0)
        t=t[w]

        s=sort(t.ellip)
        t=t[s]
        pplot, /over, t.ellip, t.pdiff, color=colors[i]

        sr2_actual = median(t.sizerat2_actual)
        add_arrval, string(sr2_actual,f='(f0.2)'), leg

    endfor

    plegend, ['admom',psfmodel+' PSF',galmodel+' galaxy'], /top, /left
    leg = textoidl('\sigma^2_{psf}/\sigma^2_{gal}')+' '+leg

    ; now reverse things
    leg=reverse(leg)
    colors_rev=reverse(colors)
	plegend,leg, /right, /top, line=0, colors=colors_rev, $
        charsize=1, spacing=1.5


    plot, [0], /nodata, $
        yrange=yrange, ystyle=3, xrange=xrange, xstyle=3, $
        xtitle='intrinsic ellip', ytitle=textoidl('\deltae/e');, $
        ;aspect=esheldon_config('gratio')

    for i=0L, nsr2-1 do begin

        t=self->getstats(psfmodel, galmodel, sr2[i])
        w=where(t.ellip gt 0)
        t=t[w]

        s=sort(t.ellip)
        t=t[s]
        pplot, /over, t.ellip, t.pdiff_rg, color=colors[i]

    endfor

    plegend, ['regauss',psfmodel+' PSF',galmodel+' galaxy'], /top, /left

    if keyword_set(dops) then begin
        endplot,/trim,/png, dpi=120
    endif else begin
        !p.charsize=1
    endelse

    !p.multi=0

end

pro regauss_sim::plot_many, psfmodel, galmodel
    sr2=self->example_sizerat2()

    for i=0L, n_elements(sr2)-1 do begin
        self->plot, psfmodel, galmodel, sr2[i]
    endfor


end


pro regauss_sim::plot, psfmodel, galmodel, sizerat2

    if not in(['gauss','exp'], galmodel) then begin
        message,'galmodel must be gauss or exp'
    endif

    st = self->getstats(psfmodel, galmodel, sizerat2)

    dir=self->dir(psfmodel)
    psfile=string($
        f='("regauss-sim-",a,"-",a,"-sizerat2-",f0.2,".eps")',$
        psfmodel,galmodel,sizerat2)
    psfile=filepath(root=dir,psfile)

	begplot,psfile,/encap,/color, $
		xsize=5.0*esheldon_config('gratio'),ysize=5

    w=where(st.ellip gt 0)
	yrange=[-0.1,0.1]
	plot, st[w].ellip, st[w].pdiff, yrange=yrange, ystyle=3, $
		xtitle='intrinsic ellip', ytitle=textoidl('\deltae/e')
	pplot,/over,st[w].ellip,st[w].pdiff+st[w].pdiff_sdev,color='blue'
	pplot,/over,st[w].ellip,st[w].pdiff-st[w].pdiff_sdev,color='blue'

    sizerat2_actual = median(st[w].sizerat2_actual)
	plegend,[galmodel,$
             textoidl('\sigma^2_{psf}/\sigma^2_{gal}')+' = '+string(sizerat2,f='(f0.2)'), $
             textoidl('actual \sigma^2_{psf}/\sigma^2_{gal}')+' = '+string(sizerat2_actual,f='(f0.2)')], $
             /right

	endplot,/trim,/png
end

pro regauss_sim::write_pbs_many, psfmodel, galmodels=galmodels
    ; write the scripts as well as some .sh files if we need
    ; to run at BNL where we can only run idl on one machine

    if n_elements(galmodels) eq 0 then galmodels=['gauss','exp']
    for i=0L, n_elements(galmodels)-1 do begin
        galmodel = galmodels[i]
        sr2=self->example_sizerat2()

        for j=0L, n_elements(sr2)-1 do begin
            self->write_pbs, psfmodel, galmodel, sr2[j], pbsfile=pbsfile
            add_arrval, pbsfile, pbsfiles
        endfor
        
    endfor

    nf = n_elements(pbsfiles)

    dir = file_dirname(pbsfiles[0])
    commands = strarr(nf)
    for i=0L, nf-1 do begin
        f = file_basename(pbsfiles[i])
        commands[i] = 'bash '+f+' &> /dev/null'
    endfor

    ; we want to use 6 processors
    nproc = 6
    plist = split_array(commands, nsplit=nproc)

    for i=0L,nproc-1 do begin
        c = *plist[i]
        fout = filepath(root=dir, 'script'+string(i, f='(i02)')+'.sh')
        print,'  Writing combined script for bach: ',fout
        openw, lun, fout, /get_lun
        for j=0L, n_elements(c)-1 do begin
            printf, lun, c[j]
        endfor
        free_lun, lun
    endfor

end
pro regauss_sim::write_pbs, psfmodel, galmodel, sizerat2, pbsfile=pbsfile

    s2str = string(sizerat2,f='(f0.2)')
    pbsdir='~/pbs/regauss-sim'
    pbsdir = filepath(root=pbsdir, psfmodel+'-psf')
    if not file_test(pbsdir) then file_mkdir, pbsdir

    pbsfile=string(f='("regauss-sim-",a,"-",a,"-",f0.2,".pbs")',psfmodel,galmodel,sizerat2)
    pbsfile=filepath(root=pbsdir, pbsfile) 
    print,'Writing pbs file: ',pbsfile

    openw, lun, pbsfile, /get_lun
    printf,lun,'#PBS -l walltime=24:00:00'
    printf,lun,'#PBS -N '+galmodel+'-'+s2str
    printf,lun,'#PBS -j oe'
    printf,lun,'#PBS -o /home/esheldon/pbs/regauss-sim/regauss-sim-'+galmodel+'-'+s2str+'.pbs.pbslog'
    printf,lun,'#PBS -m a'
    printf,lun,'#PBS -V'
    printf,lun,'#PBS -r n'
    printf,lun,'#PBS -W umask=0022'
    printf,lun,'echo Running on `hostname`'
    printf,lun
    printf,lun,'# my startup has dependencies I dont want for this'
    printf,lun,'export IDL_STARTUP=""'
    printf,lun
    printf,lun,'setup photoop'
    printf,lun,'setup idlutils'
    printf,lun,'setup sdssidl -r ~esheldon/sdssidl'
    printf,lun,'setup esidl -r ~esheldon/exports/esidl-work'
    printf,lun,'setup tree'
    printf,lun
    printf,lun,'logf='+pbsfile+'.log'
    printf,lun,'idl &> "$logf" <<EOF'
    printf,lun,'    rs=obj_new("regauss_sim")'
    printf,lun,'    rs->run, "'+psfmodel+'", "'+galmodel+'", '+s2str
    printf,lun,'EOF'

    free_lun, lun
end

function regauss_sim::example_sizerat2, less=less
    all=[0.05,0.1,0.2,0.3,1.0/3.0,0.4,0.5,1.0,2.0]
    if keyword_set(less) then begin
        return, all[1:n_elements(all)-2]
    endif else begin
        return, all
    endelse
end


function regauss_sim::random_objects_file
	dir=getenv('REGAUSSIM_DIR')
    outfile = filepath(root=dir, 'rand_field_objects.fits')
    return,outfile
end
function regauss_sim::random_objects_read
    f = self->random_objects_file()
    if not fexist(f) then message,'File not found:'+f
    print,'reading random field objects file: ',f
    str = mrdfits(f,1)
    return, str
end
pro regauss_sim::make_random_objects_list

    ; we want to choose objects from the survey
    ; with representative seeeing.
    ;
    ; randomly choose fields from the survey
    ; and read in data until we have the desired
    ; number of objects.

    outfile=self->random_objects_file()
    print,'Will write to file: ',outfile

    ; go till we get 1,000,000 objects
    ntot = 1000000
    ;ntot = 10000

    print,'reading window flist'
    window_read, flist=flist, rescore=rescore
    print,'choosing objects with score > 0.1 and psf fwhm in [1,2]'
    indx = where(flist.score GT 0.1 and flist.psf_fwhm gt 1 and flist.psf_fwhm lt 2, ct)
    flist = flist[indx]

    ; we wont' fill all of this
    ptrlist = ptrarr(ct)

    ; get random index into structure
    s=sort( randomu(seed,ct) )
    irand=0L
    n=0L

    while (n lt ntot) do begin
        print,n,ntot,f='("count: ",i0,"/",i0)'
        ; get the next random field
        i = s[irand]

        ; read this field
        tmp = sdss_read('fpObjc',flist[i].run,flist[i].camcol,flist[i].field,rerun=flist[i].rerun, $
                        taglist=['run','rerun','camcol','field','id','colc','rowc'])
        n += n_elements(tmp)
        
        ptrlist[irand] = ptr_new(tmp, /no_copy)

        irand += 1
    endwhile

    print,'combining'
    str = combine_ptrlist(ptrlist)

    print,'writing file: ',outfile
    mwrfits, str, outfile, /create
end

function regauss_sim::r0_from_mom, mom
    ; estimate of r0 from measured
    ; admom moments, where mom=ixx+iyy
    common _r0_vs_mom_block, mom_, r0_, minmom_, maxmom_

    if n_elements(mom_) eq 0 then begin
        t=self->read_r0_vs_mom(1.0, 0.0)
        mom_ = t.ixx + t.iyy
        r0_ = t.r0

        minmom_ = min(mom_)
        maxmom_ = max(mom_)
    endif

    if mom lt minmom_ or mom gt maxmom_ then begin
        message,string(mom,minmom_,maxmom_,$
                       f='("mom value ",f0.2," is out of range [",f0.2,",",f0.2,"]")')
    endif

    r0 = interpol(r0_, mom_, mom)
    return, r0
end



; determinte the relationship between exponential scale, axis ratio 
; and adaptive moments ixx,iyy,ixy and ellip
; only do it for round, just want a rough estimate for the sims

function regauss_sim::r0_vs_mom_file, aratio, theta
     outfile=string(aratio,theta,f='("admom_vs_r0_aratio",f0.2,"_theta",f0.2,".st")')
     outfile=filepath(root=getenv('ESIDL_DIR'), sub=['pro','admom','data'], outfile)
     return, outfile
end
function regauss_sim::read_r0_vs_mom, aratio, theta, hdr=hdr
     file=self->r0_vs_mom_file(aratio, theta)
     return, read_idlstruct(file, hdr=hdr)
end
function regauss_sim::aratio_vs_mom_file, r0, theta
     outfile=string(r0,theta,f='("admom_vs_aratio_r0",f0.2,"_theta",f0.2,".st")')
     outfile=filepath(root=getenv('ESIDL_DIR'), sub=['pro','admom','data'], outfile)
     return, outfile
end
function regauss_sim::read_aratio_vs_mom, r0, theta, hdr=hdr
     file=self->r0_vs_mom_file(r0, theta)
     return, read_idlstruct(file, hdr=hdr)
end


pro regauss_sim::calculate_exp_r0_from_mom, nr0=nr0


    if n_elements(nr0) eq 0 then nr0=100

    r0 = arrscl(findgen(nr0), 0.52, 3.0)
    aratio = 1.0
    theta = 0.0

    outfile=self->r0_vs_mom_file(aratio, theta)
    psfile=repstr(outfile, '.st', '.eps')

    ixx=fltarr(nr0)
    iyy=ixx
    ixy=ixx
    momerr=ixx
    rho4=ixx
    whyflag=lonarr(nr0)

    imsize=[50,50]
    ;imsize=[100,100]
    sky=0.0
    skysig=5.0
    for i=0L, nr0-1 do begin
        
        tr0=r0[i]
        make_exp, image, imsize, tr0, aratio=aratio, theta=theta, $
            cen=cen

        wguess = tr0^2
        admom, image, cen[0], cen[1], sky, skysig, wguess, $
            tixx, tiyy, tixy, tmomerr, trho4, twhyflag
        if twhyflag ne 0 then begin
            print,'whyflag not zero: ',twhyflag
        endif

        ixx[i]=tixx
        iyy[i]=tiyy
        ixy[i]=tixy
        momerr[i]=tmomerr
        rho4[i]=trho4
        whyflag[i]=twhyflag
    endfor

    if !d.name eq 'PS' then begin
        fitcolor='blue' 
    endif else begin
        ;!p.charsize=2
        fitcolor='green'
    endelse

    w=where(ixx gt 0)

    sz=ixx[w]+iyy[w]
    ;xrange=[0,1.1*max(sz)]
    ;yrange=xrange

    ;!p.multi=[0,2,0]
    ;rdis, image

    ;pplot, sz, r0[w]^2,  $
    begplot, psfile, /encap, /color
    pplot, sz, r0[w]^2/sz,  $
        xrange=xrange, yrange=yrange, psym=-8, $
        ytitle=textoidl('r_0^2/(I_{xx}+I_{yy})'), $
        xtitle=textoidl('adaptive I_{xx}+I_{yy}'), $
        aspect=1

    legend,$
        [string(aratio,f='("aratio: ",f0.2)'),$
         string(theta,f='("theta: ",f0.2)')], /bottom, /right
     endplot,/trim,/png

    ;fitpars = linfit( sz, r0[w]^2, yfit=yfit, sigma=sigma )
    ;pplot, /overplot, sz, yfit, color=fitcolor

    ;print,fitpars[0],fitpars[1],f='("r0^2 = ",f0.2," + ",f0.2," (ixx+iyy)")'

    ;legend,$
    ;    [string(aratio,f='("aratio: ",f0.2)'),$
    ;     string(theta,f='("theta: ",f0.2)'),$
    ;     textoidl('fit: r_0^2 = a + b (I_{xx}+I_{yy})'),$
    ;     string(fitpars[0],sigma[0],f='("a: ",f0.3,"+/-",f0.4)'),$
    ;     string(fitpars[1],sigma[1],f='("b: ",f0.3,"+/-",f0.4)')]


     hdr={aratio:aratio, $
         theta:theta}

     st={r0:0.0, ixx:0.0, ixy:0.0, iyy:0.0}
     st=replicate(st, nr0)
     st.r0=r0
     st.ixx=ixx
     st.ixy=ixy
     st.iyy=iyy

     print,'writing to output file: ',outfile

     write_idlstruct, st, outfile, /csv, hdr=hdr


end


pro regauss_sim::calculate_exp_aratio_from_mom, naratio=naratio

    ; ~ gaussian FHWM 2''
    ;r0=0.63
    r0=2.0

    if n_elements(naratio) eq 0 then naratio=100

    aratio = arrscl(findgen(naratio), 0.3, 1.0)
    theta = 0.0

    outfile=self->aratio_vs_mom_file(r0, theta)
    psfile=repstr(outfile, '.st', '.eps')
    print,outfile,psfile

    ixx=fltarr(naratio)
    iyy=ixx
    ixy=ixx
    momerr=ixx
    rho4=ixx
    whyflag=lonarr(naratio)

    ;imsize=[50,50]
    imsize=[100,100]
    sky=0.0
    skysig=5.0
    for i=0L, naratio-1 do begin
        
        taratio = aratio[i]
        make_exp, image, imsize, r0, aratio=taratio, theta=theta, $
            cen=cen

        wguess = r0^2
        admom, image, cen[0], cen[1], sky, skysig, wguess, $
            tixx, tiyy, tixy, tmomerr, trho4, twhyflag
        if twhyflag ne 0 then begin
            print,'whyflag not zero: ',twhyflag
        endif

        ixx[i]=tixx
        iyy[i]=tiyy
        ixy[i]=tixy
        momerr[i]=tmomerr
        rho4[i]=trho4
        whyflag[i]=twhyflag
    endfor

    w=where(ixx gt 0, nkeep)

    T = ixx[w] + iyy[w]
    e1 = (ixx[w]-iyy[w])/T
    e2 = 2.0*ixy[w]/T
    ellip = sqrt(e1^2 + e2^2)

    ;begplot, psfile, /encap, /color
    if !d.name eq 'PS' then begin
        fitcolor='blue' 
    endif else begin
        ;!p.charsize=2
        fitcolor='green'
    endelse


    pplot, ellip, aratio[w],  $
        xrange=xrange, yrange=yrange, psym=-8, $
        ytitle=textoidl('axis ratio'), $
        xtitle=textoidl('adaptive ellip'), $
        aspect=1, /ynozero

    pplot, /over, ellip, (1-ellip)/(1+ellip), color=fitcolor

    legend,$
        [string(r0,f='("r0: ",f0.2)'),$
         string(theta,f='("theta: ",f0.2)')], /bottom, /right
     endplot,/trim,/png

     hdr={r0:r0, $
         theta:theta}

     st={aratio:0.0, ixx:0.0, ixy:0.0, iyy:0.0, ellip:0.0}
     st=replicate(st, nkeep)
     st.aratio=aratio[w]
     st.ixx=ixx[w]
     st.ixy=ixy[w]
     st.iyy=iyy[w]
     st.ellip=ellip

     print,'writing to output file: ',outfile

     write_idlstruct, st, outfile, /csv, hdr=hdr


     endplot

end


pro regauss_sim::test_mom2disk_exp, theta

    ; ~ gaussian FHWM 1.5''
    T = 5.0

    nellip=100

    ellip = arrscl(findgen(nellip), 0.0, 0.4)

    st = { $
        T_input: -9999.0, $
        theta_input: -9999.0, $
        e_input: -9999.0, $
        e1_input: -9999.0, $
        e2_input: -9999.0, $
        ixx: -9999.0, $
        ixy: -9999.0, $
        iyy: -9999.0, $
        rho4: -9999.0, $
        momerr: -9999.0, $
        whyflag: 0L, $
        e: -9999.0, $
        e1: -9999.0, $
        e2: -9999.0 $
    }
    st=replicate(st, nellip)

    imsize=[100,100]
    for i=0L, nellip-1 do begin
        self->regauss::ellip2mom, ellip[i], theta, T, ixx, ixy, iyy

        disk = mom2disk('exp', ixx, ixy, iyy, imsize, counts=1.0, $
                        cen=cen)

        sky=0.0
        skysig=5.5
        wguess=T/2.0
        admom, $
            disk, cen[0], cen[1], sky, skysig, wguess, $
            tixx, tiyy, tixy, tmomerr, trho4, twhyflag

        st[i].T_input = T
        st[i].theta_input = theta

        st[i].e_input = ellip[i]
        st[i].e1_input = (ixx-iyy)/T
        st[i].e2_input = 2.0*ixy/T

        st[i].ixx = tixx
        st[i].ixy = tixy
        st[i].iyy = tiyy
        st[i].momerr = tmomerr
        st[i].rho4 = trho4
        st[i].whyflag = twhyflag

        if twhyflag eq 0 then begin
            st[i].e1 = (tixx-tiyy)/T
            st[i].e2 = 2.0*tixy/T
            st[i].e = sqrt(st[i].e1^2 + st[i].e2^2)
        endif else begin
            print,'failure: ',twhyflag
        endelse

    endfor

    w=where(st.whyflag eq 0)
    st=st[w]

    pplot, st.e_input, st.e, $
        xtitle='ellip input', $
        ytitle='mom2disk exp ellip', $
        psym=-8
    pplot,/over,st.e_input, st.e_input, line=2
end





pro regauss_sim__define

  struct = {$
             regauss_sim, $
             inherits regauss $
           }

end 
