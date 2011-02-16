
;+
; NAME:
;   goods_sensitivity__define
;     
;
; PURPOSE:
;   An IDL class to calculate the sensitivity from degraded goods images.  The
;   actual degrading process is performed in goods_addnoise__define
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
;-
function goods_sensitivity::init
    return, 1
end


;
; Files and directories
;
function goods_sensitivity::dir, subdir=subdir, makedir=makedir
    dir = esheldon_config('des_sensitivity_dir')
    dir = path_join(dir, 'sensitivity.new')
    if n_elements(subdir) ne 0 then begin
        dir = path_join(dir, subdir=subdir)
    endif

    if keyword_set(makedir) then begin
        if not fexist(dir) then file_mkdir, dir
    endif
    return, dir
end
function goods_sensitivity::file, subdir=subdir, prefix=prefix, suffix=suffix, ext=ext, makedir=makedir
    dir = self->dir(subdir=subdir, makedir=makedir)
    if n_elements(ext) eq 0 then ext='.fits'

    file_el = 'convolved'

    if n_elements(prefix) ne 0 then file_el = [prefix,file_el]
    if n_elements(suffix) ne 0 then file_el = [file_el,suffix]
    file = strjoin(file_el,'-')+ext

    file = filepath(root=dir, file)

    return, file
end

function goods_sensitivity::read_catalogs, type=type, seeing=seeing, admom=admom

	; default is read admom
	if n_elements(admom) eq 0 then admom=1

	if keyword_set(admom) then begin
		admomstr = 'admom*' 
	endif else begin
		admomstr=''
	endelse
	if n_elements(seeing) eq 0 then seeing=0.0

	if seeing gt 0.0 then begin
		seestr = ntostr(seeing, 4, /round)+'_'
		if n_elements(type) eq 0 then typestr = '' else typestr=type+'_'
	endif else begin
		seestr = ''
	endelse


    dir = esheldon_config('des_sensitivity_dir')
	catdir = path_join(dir, 'cat')

	if seeing eq 0.0 then begin
		pattern = path_join(catdir, "h_goods*"+admomstr+".fit")
	endif else begin
		pattern = path_join(catdir, typestr+seestr+"h_goods*"+admomstr+".fit")
	endelse

	;pattern = path_join(catdir, "h_goods*admom*.fit")
	print,'Searching for pattern: ',pattern
	files = findfile(pattern)
      
	print,'Reading files: ',files
	struct = mrdfits_multi(files)

	return, struct

end
;
; Calculate lensing weights
;

function goods_sensitivity::weight, shapenoise, sigma
    return, 1.0/(shapenoise^2 + sigma^2)
end
function goods_sensitivity::weight_of_zl, shapenoise, sigma, zl, zs
    weight = self->weight(shapenoise, sigma)
    siginv = sigmacritinv(zl, zs)
    weight = weight*siginv^2
    return, weight
end

;
; Convert between effective density and sensitivity
;
function goods_sensitivity::sens2neff, shapenoise, sens, area
    return, (shapenoise*0.5/sens)^2/area
end
function goods_sensitivity::neff2sens, shapenoise, neff, area
    return, shapenoise*0.5/sqrt( neff*area )
end


;
; Calculate the sensitivity from the shape distribution
;
function goods_sensitivity::calc_sens_neff, e1, e2, weights
    common _htsc_common, shapenoise, area, survey_area, nexp, minmag, maxmag, $
        nzs_sample, arcperpix

    wmom, e1, blah, we1mean, we1sig, we1err, weights=weights
    wmom, e2, blah, we2mean, we2sig, we2err, weights=weights

    sens = 0.5*mean( [we1err, we2err] )
    neff = self->sens2neff(shapenoise, sens, area)

    return, {sens:sens, neff:neff}
end

;
; Calculate the sensitivity and neff in bins of some size
;
function goods_sensitivity::calc_sens_neff_fwhm, e1,e2,weights,fwhm,fwhmbin

    common _htsc_common, shapenoise, area, survey_area, nexp, minmag, maxmag, $
        nzs_sample, arcperpix

    if n_elements(maxmag) eq 0 then maxmag = 30.0

    weights_scaled = weights/max(weights)
    bs = binner(fwhm, weights_scaled, binsize=fwhmbin, rev=rev)

    ; plot weight vs size
    xt = 'adaptive fwhm'
    yt = 'weight'
    ;xrange = [minmag, maxmag]
    pplot, bs.xmean, bs.ymean, yerr=bs.yerr, psym=8, aspect=1, $
        yrange = [0.0, 1.2], xrange=xrange, /xsty, xtit=xt, ytit=yt
    pplot, $
        [bs.xmean, max(bs.xmean)+fwhmbin/2.],  $
        [bs.hist/float(max(bs.hist))/2.,0.0],psym=10, $
        /overplot

    legend,'IAB < '+ntostr(maxmag,f='(I0.2)'), box=0

    key=prompt_kbrd()

    nbin = n_elements(bs.xmean)
    msens = replicate(9999d, nbin)
    mneff = replicate(0d, nbin)
    for i=0L, nbin-1 do begin
        if rev[i] ne rev[i+1] then begin
            w=rev[ rev[i]:rev[i+1]-1 ]
            if n_elements(w) gt 1 then begin

                tst = self->calc_sens_neff(e1[w], e2[w], weights[w])

                ; This is shear sensitivity for this field.
                msens[i] = tst.sens
                mneff[i] = tst.neff
            endif 
        endif
    endfor

    retstruct = { $
        fwhm: bs.xmean, $
        sweight: bs.ymean ,$
        sweight_err: bs.yerr, $
        shist: bs.hist, $
        sdensity:bs.hist/area, $
        ssens: msens, $
        sneff: mneff $
    }

    return, retstruct

end



;
; Calculate the sensitivity and neff in bins of magnitude
;
function goods_sensitivity::calc_sens_neff_mags, e1, e2, weights, mags, magbin

    common _htsc_common, shapenoise, area, survey_area, nexp, minmag, maxmag, $
        nzs_sample, arcperpix

    if n_elements(maxmag) eq 0 then maxmag = 30.0

    weights_scaled = weights/max(weights)
    bs = binner(mags, weights_scaled, binsize=magbin, rev=rev)

    ; plot weight vs mag
    xt = 'IAB'
    yt = 'weight'
    xrange = [minmag, maxmag]
    pplot, bs.xmean, bs.ymean, yerr=bs.yerr, psym=8, aspect=1, $
        yrange = [0.0, 1.2], xrange=xrange, /xsty, xtit=xt, ytit=yt
    pplot, $
        [bs.xmean, max(bs.xmean)+magbin/2.],  $
        [bs.hist/float(max(bs.hist))/2.,0.0],psym=10, $
        /overplot

    legend,'IAB < '+ntostr(maxmag,f='(I0.2)'), box=0

    key=prompt_kbrd()

    nbin = n_elements(bs.xmean)
    mag_auto = bs.xmean
    msens = replicate(9999d, nbin)
    mneff = replicate(0d, nbin)
    for i=0L, nbin-1 do begin
        if rev[i] ne rev[i+1] then begin
            w=rev[ rev[i]:rev[i+1]-1 ]
            if n_elements(w) gt 1 then begin

                tst = self->calc_sens_neff(e1[w], e2[w], weights[w])

                ; This is shear sensitivity for this field.
                msens[i] = tst.sens
                mneff[i] = tst.neff
            endif 
        endif
    endfor

    retstruct = { $
        mag_auto: bs.xmean, $
        mweight: bs.ymean ,$
        mweight_err: bs.yerr, $
        mhist: bs.hist, $
        mdensity:bs.hist/area, $
        msens: msens, $
        mneff: mneff $
    }

    return, retstruct

end


;
; Generate random redshifts from an input redshift distribution
;

function goods_sensitivity::generate_random_zsource, zs, pofzs, nrand
    genrand, pofzs, zs, nrand, randzs, method=1, /double
    return, randzs
end


;
; Generate random redshifts for galaxies from their magnitudes and P(z|m)
;

function goods_sensitivity::generate_random_z_from_mags, mags, pzm, magbin

    common _htsc_common, shapenoise, area, survey_area, nexp, minmag, maxmag, $
        nzs_sample, arcperpix

    nmags = n_elements(mags)
    npzm = n_elements(pzm)
    nmb = n_elements(magbin)
    if (nmags eq 0) or (npzm eq 0) or (nmb eq 0) then begin
        message,'usage: rz=obj->generate_random_z_from_mags(mag, pzm, magbin)'
    endif

    nobj = n_elements(mags)
    randz = dblarr(nobj)

    ; Bin into magnitude bins
    bs = binner(mags, min=minmag, max=maxmag, binsize=magbin, rev=rev)

    ; Now interpolate Huan's p(z|m) to each magnitude bin 
    ist = self->huan_interp(pzm, bs.xcenter)

    nbin = n_elements(bs.xcenter)
    for i=0L, nbin-1 do begin
        if rev[i] ne rev[i+1] then begin
            w=rev[ rev[i]:rev[i+1]-1 ]
            nw=n_elements(w)
            ; Generate some redshifts for this mag bin
            randz[w] = $
                self->generate_random_zsource(ist[i].z, ist[i].dndzda, nw)

            ;plothist, z[w], bin=0.1, /norm
            ;pplot, $
            ;    ist[i].z, ist[i].dndzda/qgauss(ist[i].dndzda, ist[i].z, 100), $
            ;    /over, color='red'
            ;key=prompt_kbrd('hit a key')
        endif
    endfor

    ;plothist, randz,bin=0.1, /norm, ytitle='P(zrand)',xtit='zrand'
    ;key=prompt_kbrd('hit a key')
    return, randz
    
end

;
; Calculate the effective number density as a function of zlens.  This 
; includes an explicit weighting with sigma_crit; i.e. assuming we are
; measuring \Delta\Sigma
;   outputs zl, neff(zl)
;
function goods_sensitivity::calc_neff_of_zl, zl, e1, e2, sigma, mags, pzm

    common _htsc_common, shapenoise, area, survey_area, nexp, minmag, maxmag, $
        nzs_sample, arcperpix

    ; Generate some random source redshifts given the input mags.
    magbin = 0.25

    ; Now, given these random redshifts and the input zl we can calculate
    ; the weight and effective number density

    nzl = n_elements(zl)
    sens = dblarr(nzl)
    neff = dblarr(nzl)

    for i=0L, nzl-1 do begin

        ssens = dblarr(nzs_sample)
        sneff = dblarr(nzs_sample)
        for j=0L, nzs_sample-1 do begin
            randz = self->generate_random_z_from_mags(mags, pzm, magbin)
            weights = self->weight_of_zl(shapenoise, sigma, zl[i], randz)
            stst = self->calc_sens_neff(e1, e2, weights)
            ssens[j] = stst.sens
            sneff[j] = stst.neff 
        endfor

        sens[i] = mean(ssens) 
        neff[i] = mean(sneff)

        ;tst = self->calc_sens_neff(e1, e2, weights)
        ;sens[i] = tst.sens
        ;neff[i] = tst.neff

        if i eq 0 then begin
            plothist, randz,bin=0.1, /norm, ytitle='P(zrand)',xtit='zrand',$
                title='i=0 random zs example plot', aspect=!gratio
            key=prompt_kbrd('hit a key')
        endif
        print,'.',f='(a,$)'
    endfor
    print

    return, {zl: zl, zlsens:sens, zlneff:neff}
end

;
; Calculate effective number density as a function of source redshift from
; Huan's P(z|m) data.  This is for a cosmic shear type measurement with no
; lens redshift weighting
;

function goods_sensitivity::calc_neff_of_zs, nefftot, mag_bins, mag_neff, pzm

    w=where(mag_bins gt 0)

    ; Interpolate Huan's p(z|m) to these magnitudes
    ; Then loop over redshift and calculate effective density in each
    ; redshift bin by integrating across mag axis
    ist = self->huan_interp(pzm, mag_bins[w])

    ; Now integrate over mag in each z bin 
    zs = pzm[0].z
    nzs = n_elements(zs)
    neff = dblarr(nzs)
    for j=0L, nzs-1 do begin
        ; This is a vector with an element for each magnitude
        neff[j] = qgauss( mag_neff[w]*ist.dndzda[j], mag_bins[w], 100 )
    endfor

    ; re-normalize
    neff = nefftot*neff/qgauss(neff, zs, 100)

    pplot, zs, neff, xtitle=textoidl('Z_s'), ytitle=textoidl('n_{eff}'),$
        aspect=!gratio, psym=-8
    key = prompt_kbrd('hit a key')

    return, {zs:zs, zsneff:neff}
end




;
;
; Extract Galaxies from the sextractor catalog
;
;

function goods_sensitivity::get_gals, $
		seeing, struct, rcut, minmag, maxmag, $
		wgood, wstar, wgal, r, med_seeing, $
		fwhm_range=fwhm_range, $
		use_cursor=use_cursor, $
		arcperpix=arcperpix

	if n_elements(arcperpix) eq 0 then arcperpix=0.03

    ; which aperture should we use
    ; apertures in pixels: [5, 7, 10, 15, 20]
    ; apertures in arcsec: 0.206 per pixel: [1.03, 1.442, 2.06, 3.09, 4.12]
    
    ; cuts
    badflags = 2L^0 + 2L^2 + 2L^3 + 2L^4 + 2L^5 + 2L^6 + 2L^7
    max_fwhm = 3.0

    seeing_width = 0.05
    min_fwhm = seeing - seeing_width

    if !d.name eq 'X' then begin 
        charsize=2
        star_color = c2i('green')
        bad_color = c2i('red')
    endif else begin 
        charsize=1
        star_color = c2i('blue')
        bad_color = c2i('red')
    endelse 

    xm = [20,6]
    ym = [8,4]

    pxm = [10,3]
    pym = [4,2]

    !x.margin = xm
    !y.margin = ym

    mag = struct.mag_auto

    yt = 'FWHM (arcsec)'
    xt = 'i_AB (mag_auto)'
    mrange = [15, 30]
    frange = [0, 3]

    size2 = struct.ixx + struct.iyy
    afwhm = sqrt(size2/2.)*2.35*arcperpix

    aplot, 1, mag, afwhm, psym=3, /center, xmargin=xm, ymargin=ym, $
        xrange=mrange, yrange=frange, xtitle=xt, ytitle=yt

    w=where(struct.ixx GT 0)
    w2 = where(afwhm[w] LT 2)
    help,w,w2


    ; Cuts:
    ; flags we care about: skipping the deblend, 2L^1

    print
    print,'Checking cuts'
    help,where((struct.flags AND badflags) NE 0) 
    help,where(afwhm LT min_fwhm)
    help,where(afwhm GT max_fwhm)
    help,where(mag LT minmag)
    help,where(mag GT maxmag)
    help,where(struct.whyflag NE 0)
    print

    OBJECT2_AMOMENT_FAINT =  '200000'X ; too faint for adaptive moments 
    OBJECT2_AMOMENT_UNWEIGHTED = '200000'X ; failed so tried unweighted mom
    OBJECT2_AMOMENT_SHIFT =  '400000'X ; centre moved too far while
    ;determining adaptive moments 
    OBJECT2_AMOMENT_MAXITER = '800000'X ; Too many iterations while

    help,where((struct.whyflag AND OBJECT2_AMOMENT_FAINT) NE 0)
    help,where((struct.whyflag AND OBJECT2_AMOMENT_UNWEIGHTED) NE 0)
    help,where((struct.whyflag AND OBJECT2_AMOMENT_SHIFT) NE 0)
    help,where((struct.whyflag AND OBJECT2_AMOMENT_MAXITER) NE 0)

    wgood=where( $
        (struct.flags and badflags) eq 0 and $
        afwhm gt min_fwhm and $
        afwhm lt max_fwhm and $
        struct.rho4 gt 0 and $
        mag gt minmag and $
        mag lt maxmag and $
        struct.whyflag eq 0, nw, $
        comp=bad, ncomp=nbad)

    help,struct,wgood

    if nw ne 0 then begin 
        oplot, mag[bad], afwhm[bad], $
            psym=8, color=bad_color

    ENDIF 

    ; stars
    ; Really cut for stars
    star_maxmag = 23.0
    star_minmag = 18.5
    if n_elements(fwhm_range) ne 0 then begin 
        star_min_fwhm = fwhm_range[0]
        star_max_fwhm = fwhm_range[1]
    endif else if keyword_set(use_cursor) then begin 
        box_cursor, x0, y0, nx, ny
        star_min_fwhm = y0
        star_max_fwhm = y0+ny-1
    endif else begin 
        star_max_fwhm = seeing + seeing_width
        star_min_fwhm = min_fwhm
    endelse 

    wstar=where(afwhm[wgood] LT star_max_fwhm AND $
        afwhm[wgood] GT star_min_fwhm AND $
        mag[wgood] LT star_maxmag AND $
        mag[wgood] GT star_minmag, nstar)


    wstar=wgood[wstar]
    help,wstar

    plot_box,star_minmag,star_maxmag,star_min_fwhm,star_max_fwhm,$
		color=c2i('red')

    key=prompt_kbrd()

    ;; just plot the good ones
    aplot, 1, mag[wgood], afwhm[wgood], psym=3, /center, $
        xmargin=xm, ymargin=ym, $
        xrange=mrange, yrange=frange, xtitle=xt, ytitle=yt

    plot_box,star_minmag,star_maxmag,star_min_fwhm,star_max_fwhm,$
		color=c2i('red')


    key=prompt_kbrd()

    ;; ellipticity
    ee = sqrt( struct.e1^2 + struct.e2^2 )
    e1 = struct.e1
    e2 = struct.e2
    momerr = struct.ellip_err

    ;; mean stellar size, to get rough value for rsmear
    psf_size2 = median(size2[wstar])
    med_seeing = median(afwhm[wstar])
    print
    print,'Median seeing: '+ntostr(med_seeing)
    print,'Median fwhm from cat: '+$
        ntostr(median(struct[wstar].fwhm_image)*arcperpix)

    psf_r4val = median(4./struct[wstar].rho4 - 1.)

    r = (psf_size2/size2)*psf_r4val/(4./struct.rho4 - 1.)
    wgal = where(r[wgood] LT rcut AND r[wgood] GT 0, ngal)
    wgal = wgood[wgal]
    help,wgal

    !p.multi=[0,2,2]

    !x.margin = pxm
    !y.margin = pym

    bin = 0.02
    plothist,ee[wgood],bin=bin,xtitle='ellipticity', peak=1, charsize=charsize
    plothist,ee[wstar],bin=bin,/overplot,color=star_color, peak=1
    plothist,ee[wgal],bin=bin,/overplot,color=c2i('red'), peak=1
    legend,['All','Star','Gal'],lin=[0,0,0], $
        color=[!p.color,star_color,c2i('red')],/right,box=0

    plothist,e1[wgood],bin=bin,xtitle='e!D1!N', peak=1, charsize=charsize
    plothist,e1[wstar],bin=bin,/overplot,color=star_color, peak=1
    plothist,e1[wgal],bin=bin,/overplot,color=c2i('red'), peak=1

    plothist,e2[wgood],bin=bin,xtitle='e!D2!N', peak=1, charsize=charsize
    plothist,e2[wstar],bin=bin,/overplot,color=star_color, peak=1
    plothist,e2[wgal],bin=bin,/overplot,color=c2i('red'), peak=1

    fmax = 3.0
    plothist,afwhm[wgood],bin=bin,xtitle='adaptive FWHM [arcsec]', $
        peak=1, charsize=charsize, xrange=[0.0, 2.0], /xstyle
    plothist,afwhm[wstar],bin=bin, $
        /overplot,color=star_color, peak=1
    plothist,afwhm[wgal],bin=bin, $
        /overplot,color=c2i('red'), peak=1

    !p.multi=0

    !x.margin = xm
    !y.margin = ym

    key=prompt_kbrd()

    aplot, 1, mag[wgood], r[wgood], psym=3, xmargin=xm, ymargin=ym, $
        xrange=mrange, yrange=[0,2], charsize=charsize, $
        xtitle=xt, ytitle='rsmear', /center

    oplot, mag[wstar], r[wstar], psym=8, color=star_color, $
        symsize=0.25
    oplot, [0, 50], [1,1]
    oplot, [0, 50], [rcut, rcut],color=c2i('red')

    print,mean(ee[wstar]),mean(ee[wgal])
    print,mean(e1[wstar]),mean(e1[wgal])
    print,mean(e2[wstar]),mean(e2[wgal])
    print,median(r[wstar]),median(r[wgal])

    !x.margin = pxm
    !y.margin = pym

    newstruct = {$
        oldstruct:struct, $
        seeing:seeing, $
        rcut:rcut, $
        minmag:minmag, $
        maxmag:maxmag, $
        wgood:wgood, $
        wstar:wstar, $
        wgal:wgal, $
        r:r, $
		fwhm:afwhm,$
        med_seeing: med_seeing $
    }
 
    return, newstruct
end 

;
; The ratio of rsum_m/rsum_n gives the relative sensitivity sens_m/sens_n
; note the inversion.  Need to write this up.
;
function goods_sensitivity::calc_relative_wsum, nexp, shapeerr, shapenoise  
    sum = total( 1.0/(1.0 + (1.0/nexp)*(shapeerr/shapenoise)^2) )
    return, sum
end

function goods_sensitivity::calc_multi_factors, nexp, ellip_err, shapenoise

    n_nexp = n_elements(nexp)
    st = { $
        multi_sens_fac:dblarr(n_nexp), $
        multi_neff_fac:dblarr(n_nexp) $
    }

    sum1 = self->calc_relative_wsum(1.0, ellip_err, shapenoise)
    for i=0L, n_nexp-1 do begin
        sumn = self->calc_relative_wsum(nexp[i], ellip_err, shapenoise)
        sensfac = sum1/sumn
        st.multi_sens_fac[i] = sum1/sumn
        st.multi_neff_fac[i] = 1.0/sensfac^2
    endfor

    return, st
end



function goods_sensitivity::calc_weights_and_dens, struct

    common _htsc_common, shapenoise, area, survey_area, nexp, minmag, maxmag, $
        nzs_sample, arcperpix

    ; Just galaxies
    wgal = struct.wgal

    wgal2 = where( struct.oldstruct[wgal].mag_auto lt maxmag, ngal)
    wgal = wgal[wgal2]
    mags = double( struct.oldstruct[wgal].mag_auto )

	fwhm = double( struct.oldstruct[wgal].ixx + struct.oldstruct[wgal].iyy )
	fwhm = sqrt( fwhm/2.0 )*2.35*arcperpix

    r = double( struct.r[wgal] )
    rfac = 1.0/(1.0 - r)
    ellip_err = double( struct.oldstruct[wgal].ellip_err*rfac )
    e1 = double( struct.oldstruct[wgal].e1*rfac )
    e2 = double( struct.oldstruct[wgal].e2*rfac )

    weights = self->weight(shapenoise, ellip_err)

    ; raw galaxy density
    raw_density = ngal/area

    print,'Raw density: ',raw_density

    ;
    ; Overall density
    ;

    tst = self->calc_sens_neff(e1, e2, weights)
    shear_sens = tst.sens
    neff = tst.neff

    ; weighted with shapenoise only
    shear_sens_shonly = $
        0.5*sqrt(  total( weights^2*shapenoise^2 )/total( weights )^2  )
    ;effective_density_shonly = (shapenoise*0.5/shear_sens_shonly)^2/area
    neff_shonly = $
        self->sens2neff(shapenoise, shear_sens_shonly, area)

    print,'Shear sensitivity = '+ntostr(shear_sens)
    print,'Just using shapenoise: '+ntostr(shear_sens_shonly)
    print,'Effective density for sh=',shapenoise,': '+ntostr(neff)

    ;
    ; Now the same but as a function of magnitude
    ; Mean weight as a function of magnitude. Note,e1,e2,weights are all
    ; already from the wgal index list
    ;

    print
    print,'Calculating neff(m)'
    magbin = 0.25d
    magstruct = self->calc_sens_neff_mags(e1, e2, weights, mags, magbin)

    print
    print,'Calculating neff(fwhm)'
	fwhm_bin = 0.1d
	fwhmstruct = self->calc_sens_neff_fwhm(e1, e2, weights, fwhm, fwhm_bin)

    ;
    ; Calculate neff as a function of source redshift, with no explicit
    ; lens redshift weighting, using Huan's P(z|m)
    ;

    print
    print,'Calculating neff(zs)'

    pzm = self->huan_dndz_read()
    zsstruct = self->calc_neff_of_zs($
        neff, magstruct.mag_auto, magstruct.mneff, pzm)

    ;
    ;
    ; And now as a function of lens redshift using Huan's P(z|m)
    ;
    ;

    print
    print,'Calculating neff(zl)'

    nzl = 100
    minzl = 0.1d
    maxzl = 1.2d
    zl = arrscl( dindgen(nzl), minzl, maxzl )
    zlstruct = self->calc_neff_of_zl(zl, e1, e2, ellip_err, mags, pzm)


    ;
    ; more correct way to estimate error with new exposures.  This calculates
    ; the expected relative error for n equivalent images compared to this
    ; one image
    ;

    mst = self->calc_multi_factors(nexp, ellip_err, shapenoise)

    ; Adding new exposures.  I think this old way is wrong
    nexp3 = 3
  
    shear_sens3_squared = $
        shear_sens^2/nexp3 + shear_sens_shonly^2*(1.0 - 1.0/nexp3)
    shear_sens3 = sqrt( shear_sens3_squared )
    neff3 = (shapenoise*0.5/shear_sens3)^2/area


    ; neffective stuff for all
    nstruct = {$
        seeing: double(struct.seeing), $
        med_seeing: double(struct.med_seeing), $
        shapenoise:shapenoise, $
        area:area, $
        $
        raw_density: raw_density, $
        $
        shear_sens: shear_sens, $
        neff:neff, $
        $
        shear_sens_shonly: shear_sens_shonly, $
        neff_shonly: neff_shonly, $
        $
        shear_sens3: shear_sens3, $
        neff3: neff3, $
        $
        multi_nexp: nexp, $
        multi_sens_fac: mst.multi_sens_fac, $
        multi_neff_fac: mst.multi_neff_fac, $
        $
        mag_auto: magstruct.mag_auto, $
        mweight: magstruct.mweight ,$
        mweight_err: magstruct.mweight_err, $
        mhist: magstruct.mhist, $
        mdensity: magstruct.mdensity, $
        mshear_sens:magstruct.msens, $
        mneff:magstruct.mneff, $
        $
		fwhm: fwhmstruct.fwhm, $
		sweight: fwhmstruct.sweight, $
		sweight_err: fwhmstruct.sweight_err, $
		shist: fwhmstruct.shist, $
		sdensity: fwhmstruct.sdensity, $
		ssens: fwhmstruct.ssens, $
		sneff: fwhmstruct.sneff, $
		$
        zs: double(zsstruct.zs), $
        zsneff: zsstruct.zsneff, $
        $
        zl: zlstruct.zl, $
        zlsens: zlstruct.zlsens, $
        zlneff: zlstruct.zlneff  $
    }

    return, nstruct
end

function goods_sensitivity::calc, seeing=seeing, $
        type=type, $
        minmag=minmag_in, maxmag=maxmag_in, $
        struct=struct, $
        wgood=wgood, $
        wstar=wstar, $
        wgal=wgal, $
        r=r, $
        fwhm_range=fwhm_range,$
        use_cursor=use_cursor

    common _htsc_common, shapenoise, area, survey_area, nexp, minmag, maxmag, $
        nzs_sample, arcperpix

	if n_elements(seeing) eq 0 then seeing=0.0

    minmag = 17d
    maxmag = 30d
	arcperpix = 0.03

    if n_elements(minmag_in) ne 0 then minmag=minmag_in
    if n_elements(maxmag_in) ne 0 then maxmag=maxmag_in


    ; How many realizations of source redshift distibution?
    nzs_sample = 10

    shapenoise=0.32d
    area = 61.2958d
    ; This is for dark energy survey
    survey_area = 5000.0d
    ; Different equivalent exposures
    nexp = double( [1.0,1.5,2.0,2.5,3.0] )

    ; for extracting galaxies
    magstr = ntostr(maxmag,4)
    rcut = 0.9d

	if seeing gt 0.0 then begin
		seestr = ntostr(seeing, 4, /round)+'_'
		if n_elements(type) eq 0 then typestr = '' else typestr=type+'_'
	endif else begin
		seestr = ''
	endelse

    dir = esheldon_config('des_sensitivity_dir')

    if n_elements(struct) eq 0 then begin 
        catdir = path_join(dir, 'cat')
		if seeing eq 0.0 then begin
			pattern = path_join(catdir, "h_goods*admom*.fit")
		endif else begin
			pattern = path_join(catdir, typestr+seestr+"h_goods*admom*.fit")
		endelse
		print,'Searching for pattern: ',pattern
        files = findfile(pattern)
      
		print,'Reading files: ',files
        struct = mrdfits_multi(files)
		if size(struct,/tname) ne 'STRUCT' then begin
			print,'Error reading files: ',files
			message,'Halting'
		endif
    endif 

    seeing_guess = sqrt(0.1d^2 + seeing^2)
    maxmag_getgal = 30.0
    newstruct = self->get_gals($
        seeing_guess, struct, rcut, minmag, maxmag_getgal,$
        fwhm_range=fwhm_range,$
        use_cursor=use_cursor, $
		arcperpix=arcperpix)

    key=prompt_kbrd()

    sens_struct = self->calc_weights_and_dens(newstruct)
    return, sens_struct

end 



pro goods_sensitivity::runcalc, maxmag=maxmag, dops=dops, type=type

    if n_elements(type) eq 0 then type='des5yr'
	; use type='des5yr' for des predictions

    if n_elements(maxmag) eq 0 then maxmag = 30d
    magsuffix = 'maxmag'+ntostr(maxmag,f='(F0.1)')

    fitsfile = self->file(prefix=type, suffix=magsuffix, ext='.fits')

;  seeing = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, $
;            0.50, 0.55, 0.60, 0.65, 0.70, 0.76, $
;            0.80, 0.85, 0.90, 0.95, 1.00, 1.05, $
;            1.10, 1.15, 1.20]

    ;; for 0.7 use 0.675
    ;; for 0.8 use 0.78
    ;; for 0.9 use 0.87

	; we include the zero added seeing, zero added noise case as a reference
	seeing = [0.0, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, $
			  0.50, 0.55, 0.60, 0.65, $
			  0.675, $
			  0.70, 0.76, $
			  0.78, $
			  0.80, 0.85, $
			  0.87, $
			  0.90, 0.95, 1.00, 1.05, $
			  1.10, 1.15, 1.20]

    nseeing = n_elements(seeing)
    for i=0l, nseeing-1 do begin 

        suffix = magsuffix+'-fwhm'+ntostr(seeing[i], f='(f0.3)')
        psfile = self->file(subdir='plots', prefix=type, suffix=suffix, $
                            ext='.ps')

        if keyword_set(dops) then begplot,psfile,/color
        ;tsens_struct = $
        ;    huan_goods_sensitivity_convolved(seeing[i], type=type)
        tsens_struct = self->calc(seeing=seeing[i], maxmag=maxmag, type=type)
        if keyword_set(dops) then endplot

        if i eq 0 then begin 
            sens_struct = replicate(tsens_struct, nseeing)
        endif else begin 
            tmp = sens_struct[i]
            struct_assign, tsens_struct, tmp
            sens_struct[i] = tmp
        endelse 

    endfor 

    print
    print,'Writing to fits file: ', fitsfile
    mwrfits, sens_struct, fitsfile, /create

end 

function goods_sensitivity::huan_dndz_dir, subdir=subdir, makedir=makedir
    dir = esheldon_config('des_sensitivity_dir')
    dir = path_join(dir, 'galaxy_distribution_model_v1.0')

    if n_elements(subdir) ne 0 then begin
        dir = path_join(dir, subdir=subdir)
    endif

    if keyword_set(makedir) then begin
        if not fexist(dir) then file_mkdir, dir
    endif

    return, dir
end
function goods_sensitivity::huan_dndz_file, subdir=subdir, suffix=suffix, ext=ext,  makedir=makedir

    if n_elements(ext) eq 0 then ext='.fits'

    dir=self->huan_dndz_dir(subdir=subdir, makedir=makedir)
    file='lfmodel-out-zdist-hdf-IAB-combined'
    if n_elements(suffix) ne 0 then file=file+'-'+suffix
    file=file+ext
    file=path_join(dir,file)
    return, file
end
function goods_sensitivity::huan_dndz_read, subdir=subdir, suffix=suffix
    file = self->huan_dndz_file(suffix=suffix)
    print,'Reading file: ',file
    st=mrdfits(file,1)
    return,st
end

; Convert huan's many files to one big FITS file
pro goods_sensitivity::huan_dndz_convert
    dir=self->huan_dndz_dir()
    outfile = self->huan_dndz_file()

    minmag = [20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0]
    maxmag = [21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0]

    minstr = ntostr(minmag, f='(f0.2)')
    maxstr = ntostr(maxmag, f='(f0.2)')

    files = path_join(dir,'lfmodel.out.zdist.hdf.IAB'+minstr+'_'+maxstr)
    stdef={z:0.0,dndzda:0.0,fraction:0.0,dndzda_cum:0.0,fraction_cum:0.0}  
    read_struct, files[0], stdef, tmpst, skipline=4

    arrval = dblarr(n_elements(tmpst))
    stdefnew = { $
        minmag:0.0, meanmag:0.0, maxmag:0.0, $
        z:arrval, $
        dndzda: arrval, $
        fraction: arrval, $
        dndzda_cum: arrval, $
        fraction_cum: arrval $
    }

    nf = n_elements(files)
    st = replicate(stdefnew, nf)

    psfile = path_join(dir,'lfmodel-out-zdist-hdf-IAB-combined.eps')
    begplot, psfile,/encap,/color
    colors=make_rainbow(nf)

    for i=0L, nf-1 do begin
        read_struct, files[i], stdef, tmpst, skipline=4

        st[i].minmag = minmag[i]
        st[i].meanmag = (minmag[i]+maxmag[i])/2
        st[i].maxmag = maxmag[i]

        st[i].z = tmpst.z
        st[i].dndzda = tmpst.dndzda
        st[i].fraction = tmpst.fraction
        st[i].dndzda_cum = tmpst.dndzda_cum
        st[i].fraction_cum = tmpst.fraction_cum


        w=where(st[i].fraction gt 0)
        if i eq 0 then begin
            pplot, st[i].z[w], st[i].fraction[w], psym=10, $
                xrange=[0,4.0], yrange=[0,0.5], aspect=!gratio, $
                xtitle='z', ytitle='F(z)'
        endif else begin
            pplot, st[i].z[w], st[i].fraction[w], psym=10, /overplot, $
                color=colors[i]
        endelse
    endfor

    mess = minstr + ' < iab < ' + maxstr
    legend, mess, line=0, color=colors, /right, charsize=1

    endplot,/trim

    print
    print,'Writing to file: ',outfile
    mwrfits, st, outfile, /create

end

; Interpolate huan's dndz in mag bins

function goods_sensitivity::huan_interp, pzm, mags
    nz = n_elements(pzm[0].z)
    nmag = n_elements(mags)

    ist = {mag:0d, z: double(pzm[0].z), dndzda:dblarr(nz)}
    ist = replicate(ist, nmag)
    ist.mag = mags

    ; interpolate each of the z bins across magnitude
    z = pzm[0].z
    for i=0L, nz-1 do begin        
        ist.dndzda[i] = interpol(pzm.dndzda[i], pzm.meanmag, mags)
    endfor

    return, ist
end










pro goods_sensitivity::doplots, type=type, maxmag=maxmag

    if n_elements(type) eq 0 then type = 'des5yr'

    if n_elements(maxmag) eq 0 then maxmag = 30d
    magsuffix = 'maxmag'+ntostr(maxmag,f='(F0.1)')

    fitsfile = self->file(prefix=type, suffix=magsuffix)
    print
    print,'Reading file: ',fitsfile
    t=mrdfits(fitsfile,1)

    nseeing = n_elements(t)

    ;
    ; neff vs zs
    ;


    suffix = magsuffix+'-'+'neff-vs-zs'
    psfile=self->file(subdir='plots', prefix=type, suffix=suffix,$
                      ext='.eps')
    begplot,psfile,/encap, /color, xsize=8.5, ysize=8.5
    !p.charsize=2

    xtitle = textoidl('Z_{source}')
    ytitle = textoidl('n_{eff}')
    colors = make_rainbow(nseeing)
    for i=0L, nseeing-1 do begin
        if i eq 0 then begin
            pplot, t[i].zs, t[i].zsneff, $
                xtitle=xtitle, ytitle=ytitle, $
                xrange=xrange, xstyle=3, aspect=!gratio, $
                yrange=yrange
        endif
        pplot, t[i].zs, t[i].zsneff, /over, color=colors[i]

    endfor

    mess = ntostr(t.seeing, f='(F0.2)')
    legend, mess, /right, charsize=0.8, line=0, color=colors

    endplot, /trim


    
    ;
    ; neff vs zl
    ;
    suffix = magsuffix+'-'+'neff-vs-zl'
    psfile=self->file(subdir='plots', prefix=type, suffix=suffix,$
                      ext='.eps')
    begplot,psfile,/encap, /color, xsize=8.5, ysize=8.5
    !p.charsize=2
    xtitle = textoidl('Z_{lens}')
    ytitle = textoidl('n_{eff}')
    colors = make_rainbow(nseeing)
    xrange = [0.1, 1.5]
    yrange = [0.2, 50]
    for i=0L, nseeing-1 do begin
        if i eq 0 then begin
            pplot, t[i].zl, t[i].zlneff, $
                xtitle=xtitle, ytitle=ytitle, $
                xrange=xrange, xstyle=3, aspect=!gratio, $
                /ylog, yrange=yrange, ystyle=3
        endif
        pplot, t[i].zl, t[i].zlneff, /over, color=colors[i]

    endfor

    mess = ntostr(t.seeing, f='(F0.2)')
    legend, mess, /right, charsize=0.8, line=0, color=colors

    endplot, /trim



    ; Plot of total sensitivity, effective density vs. seeing
    suffix = magsuffix+'-'+'neff-vs-seeing'
    psfile=self->file(subdir='plots', prefix=type, suffix=suffix,$
                      ext='.eps')
    begplot,psfile,/encap,/color

    xrange = [0.15, 1.3]
    xtitle='seeing'
    erase & multiplot, [1,2], /square, mXtitle=xtitle, $
        mxtitsize=1.5, mxTitOffset=1.0

    multi_colors = [!p.color, c2i(['blue','red','darkgreen', 'magenta'])]

    ytitle='Sensitivity'
    pplot, t.seeing, t.shear_sens, psym=8, $
        xrange=xrange, xstyle=3, ytitle=ytitle
    ; with two equivalent bands
    nn=n_elements(t[0].multi_nexp)
    for i=1,nn-1 do begin
        new_sens = t.shear_sens*t.multi_sens_fac[i]
        pplot, t.seeing, new_sens, /overplot, color=multi_colors[i]
    endfor
    mess=ntostr(t[0].multi_nexp,f='(f0.1)')
    mess[0] = 'nexp='+mess[0]
    legend, mess, $
        psym=8, color=multi_colors, /right, /bottom

    multiplot

    ytitle='Effective Density'
    yrange = [0,30]
    pplot, t.seeing, t.neff, psym=8, $
        xrange=xrange, xstyle=3, ytitle=ytitle, $
        yrange = yrange, ystyle=1

    for i=1,nn-1 do begin
        new_neff = t.neff*t.multi_neff_fac[i]
        pplot, t.seeing, new_neff, /overplot, color=multi_colors[i]
    endfor

    multiplot,/default

    endplot,/trim_bbox


    ;
    ; Relative weight vs mag
    ; plot for seeing=0.9
    ;

    n=16
    seeing=ntostr(t[n].med_seeing, f='(f0.2)')

    suffix = magsuffix+'-seeing'+seeing+'-weight-vs-mag'
    psfile=self->file(subdir='plots', prefix=type, suffix=suffix,$
                      ext='.eps')
    begplot,psfile,/encap
    xrange=[20,26]
    yrange=[0,1.1]
    xtitle='i AB'
    ytitle='Relative Weight'
    pplot, t[n].mag_auto, t[n].mweight, yerr=t[n].mweight_err, psym=8, $
        aspect=1, xrange=xrange, xstyle=3, yrange=yrange, ystyle=1, $
        xtitle=xtitle, ytitle=ytitle
    w=where(t[n].mhist gt 0)
    hist = t[n].mhist[w]*1.0/max(t[n].mhist[w])*0.5
    pplot, t[n].mag_auto[w], hist, psym=10, $
        /overplot
    legend,'seeing='+seeing,/right
    endplot,/trim


    ;
    ; Effective density vs mag
    ; plot for seeing=0.9
    ;
    suffix = magsuffix+'-seeing'+seeing+'-density-vs-mag'
    psfile=self->file(subdir='plots', prefix=type, suffix=suffix,$
                      ext='.eps')
    begplot,psfile,/encap, /color

    w=where(t[n].mhist gt 0)
    ytitle='Number Density'
    xrange=[19,26]
    yrange=[0,30]
    pplot, t[n].mag_auto[w], t[n].mdensity[w], psym=10, aspect=1, $
        xrange=xrange, xstyle=3, xtitle=xtitle
    pplot, t[n].mag_auto[w], t[n].mneff[w], $
        psym=10, /overplot, color='blue'
    legend, ['Density','Effective Density'], line=0, $
        color=[!p.color, c2i('blue')]
    legend,'seeing='+seeing,/center,/left

    endplot,/trim

    
end

; make some text files for fabian
pro goods_sensitivity::for_fabian, type=type, seeing=seeing, maxmag=maxmag

	if n_elements(type) eq 0 then type = 'des5yr'
	if n_elements(seeing) eq 0 then seeing=0.9
    if n_elements(maxmag) eq 0 then maxmag = 30d

    magsuffix = 'maxmag'+ntostr(maxmag,f='(F0.1)')

    fitsfile = self->file(prefix=type, suffix=magsuffix)
    print
    print,'Reading file: ',fitsfile
    t=mrdfits(fitsfile,1)

	indir = file_dirname(fitsfile)

	outdir = path_join(indir, 'for-fabian')
	mag_psfile = path_join(outdir, 'magdist.eps')
	fwhm_psfile = path_join(outdir, 'fwhmdist.eps')


	izero = 0
	; get closest to seeing of 0.9
	diff = abs(t.seeing-0.9)
	mdiff = min(diff,i09)


	begplot, mag_psfile, /color, /encapsulated

	magstr = 'i_{AB}'
	xtitle=textoidl(magstr)
	ytitle=textoidl('dn/d( '+magstr+' ) [arcmin^{-2}]')
	pplot, t[izero].mag_auto, t[izero].mdensity, /ylog, $
		aspect=1, ytickf='loglabels', $
		xtitle=xtitle, ytitle=ytitle
	pplot, t[i09].mag_auto, t[i09].mdensity, /overplot, $
		color='darkgreen'
	pplot, t[i09].mag_auto, t[i09].mneff, /overplot, $
		color='darkgreen', line=2

	legend, ['no seeing added', $
		'seeing='+ntostr(t[i09].seeing,f='(f0.2)'), $
		'seeing='+ntostr(t[i09].seeing,f='(f0.2)')+', weighted'], $
		color=[!p.color,c2i('darkgreen'), c2i('darkgreen')], $
		line=[0,0,2]
	endplot, /trim


	begplot, fwhm_psfile, /color, /encapsulated

	sstr = textoidl('fwhm_{AM}')
	xtitle=textoidl(sstr +' [arcsec]')
	ytitle=textoidl('dn/d( '+sstr+' ) [arcmin^{-2}]')
	pplot, t[izero].fwhm, t[izero].sdensity, /ylog, $
		aspect=1, ytickf='loglabels', $
		xtitle=xtitle, ytitle=ytitle
	pplot, t[i09].fwhm, t[i09].sdensity, /overplot, $
		color='darkgreen'
	pplot, t[i09].fwhm, t[i09].sneff, /overplot, $
		color='darkgreen', line=2

	legend, ['no seeing added', $
		'seeing='+ntostr(t[i09].seeing,f='(f0.2)'), $
		'seeing='+ntostr(t[i09].seeing,f='(f0.2)')+', weighted'], $
		color=[!p.color,c2i('darkgreen'), c2i('darkgreen')], $
		line=[0,0,2],/right
	endplot, /trim


	; output ascii files
	w=[izero, i09]
	names = ntostr(t[w].seeing, f='(f0.2)')

	magstr_def = {iab:0d, density:0d, neff:0d}
	fwhmstr_def = {fwhm:0d, density:0d, neff:0d}
	for i=0L, n_elements(names)-1 do begin
		ii = w[i]
		name = names[i]

		magstr = replicate(magstr_def, n_elements(t[ii].mag_auto) )
		fwhmstr = replicate(fwhmstr_def, n_elements(t[ii].fwhm) )

		magstr.iab = t[ii].mag_auto
		magstr.density = t[ii].mdensity
		magstr.neff = t[ii].mneff

		fwhmstr.fwhm = t[ii].fwhm
		fwhmstr.density = t[ii].sdensity
		fwhmstr.neff = t[ii].sneff

		hdr={seeing:t[ii].seeing}

		magfile = path_join(outdir, 'magdist-'+name+'.tab')
		print,'Writing file: ',magfile
		write_idlstruct, magstr, magfile, /ascii, hdr=hdr


		fwhmfile = path_join(outdir, 'fwhmdist-'+name+'.tab')
		print,'Writing file: ',fwhmfile
		write_idlstruct, fwhmstr, fwhmfile, /ascii, hdr=hdr
	endfor



	; pre-seeing and noise catalogs
	outstdef = {$
		isgal:0, $
		mag_auto:9999d, $
		fwhm:-9999d, $
		r:-9999d, $
		ellip_err:9999d, $
		weight:0d}

	galfile = path_join(outdir, 'size-weight.tab')
	galpsfile = path_join(outdir, 'size-weight.ps')
	struct = self->read_catalogs()

	seeing_none=0.0
    seeing_guess = sqrt(0.1d^2 + seeing_none^2)
	minmag=19.0
	maxmag=30.0
	rcut=0.9d
	begplot, galpsfile, /color
	gstruct = self->get_gals(seeing_guess, struct, rcut, minmag, maxmag)
	endplot

	wgal = gstruct.wgal
	;outstruct = replicate(outstdef, n_elements(wgal))
	outstruct = replicate(outstdef, n_elements(gstruct.oldstruct))

	outstruct[wgal].isgal = 1

	r = double( gstruct.r[wgal] )
    rfac = 1.0/(1.0 - r)
    ellip_err = double( gstruct.oldstruct[wgal].ellip_err*rfac )
    ;ellip_err = gstruct.oldstruct.ellip_err*rfac )
	;ellip_err[wgal] = ellip_err[wgal]*rfac

	weights = self->weight(0.32, ellip_err)

	outstruct.mag_auto = gstruct.oldstruct.mag_auto

	outstruct.fwhm = gstruct.fwhm
	w=where(outstruct.fwhm ne outstruct.fwhm, nw)
	if nw ne 0 then outstruct[w].fwhm = -9999

	outstruct.r = gstruct.r

	outstruct[wgal].ellip_err = ellip_err[wgal]
	outstruct[wgal].weight = weights

	print,'Writing to file: ',galfile
	write_idlstruct, outstruct, galfile, /ascii




	; seeing
    seeing_guess = sqrt(0.1d^2 + seeing^2)
	seeingstr = ntostr(seeing_guess, 4, /round)

	pel = ['size-weight',type,seeingstr]
	pel = strjoin(pel, '-')
	galfile = path_join(outdir, pel+'.tab')
	galpsfile = path_join(outdir, pel+'.ps')

	struct = self->read_catalogs(type=type, seeing=seeing)

	minmag=19.0
	maxmag=30.0
	rcut=0.9d
	begplot, galpsfile, /color
	gstruct = self->get_gals(seeing_guess, struct, rcut, minmag, maxmag)
	endplot

	wgal = gstruct.wgal
	outstruct = replicate(outstdef, n_elements(gstruct.oldstruct))

	outstruct[wgal].isgal = 1

	r = double( gstruct.r[wgal] )
    rfac = 1.0/(1.0 - r)
    ellip_err = double( gstruct.oldstruct[wgal].ellip_err*rfac )

	weights = self->weight(0.32, ellip_err)

	outstruct.mag_auto = gstruct.oldstruct.mag_auto

	outstruct.fwhm = gstruct.fwhm
	w=where(outstruct.fwhm ne outstruct.fwhm, nw)
	if nw ne 0 then outstruct[w].fwhm = -9999

	outstruct.r = gstruct.r

	outstruct[wgal].ellip_err = ellip_err
	outstruct[wgal].weight = weights

	print,'Writing to file: ',galfile
	write_idlstruct, outstruct, galfile, /ascii



end

pro goods_sensitivity__define
    struct = {goods_sensitivity, _goods_sens_dummy:0}
end
