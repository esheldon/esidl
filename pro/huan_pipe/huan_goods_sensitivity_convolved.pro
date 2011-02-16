function hgsc_get_gals, seeing, struct, rcut, minmag, maxmag, $
                                   wgood, wstar, wgal, r, med_seeing, $
                                   fwhm_range=fwhm_range, $
                                   use_cursor=use_cursor

    ; which aperture should we use
    ; apertures in pixels: [5, 7, 10, 15, 20]
    ; apertures in arcsec: 0.206 per pixel: [1.03, 1.442, 2.06, 3.09, 4.12]
    
    ; cuts
    badflags = 2L^0 + 2L^2 + 2L^3 + 2L^4 + 2L^5 + 2L^6 + 2L^7
    max_fwhm = 3.0

    seeing_width = 0.05
    min_fwhm = seeing - seeing_width

    arcperpix = 0.03              ; !!!

    if !d.name eq 'X' then begin 
        charsize=2
        star_color = !green
        bad_color = !red
    endif else begin 
        charsize=1
        star_color = !blue
        bad_color = !red
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

    plot_box,star_minmag,star_maxmag,star_min_fwhm,star_max_fwhm,color=!red

    key=prompt_kbrd()

    ;; just plot the good ones
    aplot, 1, mag[wgood], afwhm[wgood], psym=3, /center, $
        xmargin=xm, ymargin=ym, $
        xrange=mrange, yrange=frange, xtitle=xt, ytitle=yt

    plot_box,star_minmag,star_maxmag,star_min_fwhm,star_max_fwhm,color=!red

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
    plothist,ee[wgal],bin=bin,/overplot,color=!red, peak=1
    legend,['All','Star','Gal'],lin=[0,0,0], $
        color=[!p.color,star_color,!red],/right,box=0

    plothist,e1[wgood],bin=bin,xtitle='e!D1!N', peak=1, charsize=charsize
    plothist,e1[wstar],bin=bin,/overplot,color=star_color, peak=1
    plothist,e1[wgal],bin=bin,/overplot,color=!red, peak=1

    plothist,e2[wgood],bin=bin,xtitle='e!D2!N', peak=1, charsize=charsize
    plothist,e2[wstar],bin=bin,/overplot,color=star_color, peak=1
    plothist,e2[wgal],bin=bin,/overplot,color=!red, peak=1

    fmax = 3.0
    plothist,afwhm[wgood],bin=bin,xtitle='adaptive FWHM [arcsec]', $
        peak=1, charsize=charsize, xrange=[0.0, 2.0], /xstyle
    plothist,afwhm[wstar],bin=bin, $
        /overplot,color=star_color, peak=1
    plothist,afwhm[wgal],bin=bin, $
        /overplot,color=!red, peak=1

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
    oplot, [0, 50], [rcut, rcut],color=!red

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
        med_seeing: med_seeing $
    }
 
    return, newstruct
end 

; The ratio of rsum_m/rsum_n gives the relative sensitivity sens_m/sens_n
; note the inversion.  Need to write this up.
function hgsc_calc_relative_wsum, nexp, shapeerr, shapenoise  
    sum = total( 1.0/(1.0 + (1.0/nexp)*(shapeerr/shapenoise)^2) )
    return, sum
end

function hgsc_calc_multi_factors, nexp, ellip_err, shapenoise

    n_nexp = n_elements(nexp)
    st = { $
        multi_sens_fac:fltarr(n_nexp), $
        multi_effdens_fac:fltarr(n_nexp) $
    }

    sum1 = hgsc_calc_relative_wsum(1.0, ellip_err, shapenoise)
    for i=0L, n_nexp-1 do begin
        sumn = hgsc_calc_relative_wsum(nexp[i], ellip_err, shapenoise)
        sensfac = sum1/sumn
        st.multi_sens_fac[i] = sum1/sumn
        st.multi_effdens_fac[i] = 1.0/sensfac^2
    endfor

    return, st
end

function hgsc_calc_weights_and_dens, struct

    common _htsc_common, shapenoise, area, survey_area, nexp

    ; Just galaxies
    wgal = struct.wgal
    ngal = n_elements(wgal)
    r = struct.r[wgal]
    rfac = 1.0/(1.0 - r)
    ellip_err = struct.oldstruct[wgal].ellip_err*rfac
    e1 = struct.oldstruct[wgal].e1*rfac
    e2 = struct.oldstruct[wgal].e2*rfac
    weights = 1./( shapenoise^2 + ellip_err^2 )

    ; raw galaxy density
    raw_density = ngal/area

    ; Now calculate the effective density
    weights_scaled = weights/max(weights)

    ngal_weighted = round( total(weights_scaled) )
    weighted_density = ngal_weighted/area

    print,'Raw density: ',raw_density
    print,'Weighted density: ',weighted_density

    ;
    ; Overall density
    ;

    wmom, e1, blah, we1mean, we1sig, we1err, weights=weights
    wmom, e2, blah, we2mean, we2sig, we2err, weights=weights

    ; This is shear sensitivity for this field.
    shear_sens = 0.5*mean([we1err, we2err])
    effective_density = (shapenoise*0.5/shear_sens)^2/area

    ; Shear sensitivity using only shape noise: Just do shapenoise*0.5/sqrt(n)
    ;shear_sens_shonly = 0.5*shapenoise/sqrt(ngal)
    ;effective_density_shonly = (shapenoise*0.5/shear_sens_shonly)^2/area

    ; weighted with shapenoise only
    shear_sens_shonly = $
        0.5*sqrt(  total( weights^2*shapenoise^2 )/total( weights )^2  )
    effective_density_shonly = (shapenoise*0.5/shear_sens_shonly)^2/area

    print,'Shear sensitivity = '+ntostr(shear_sens)
    print,'Just using shapenoise: '+ntostr(shear_sens_shonly)
    print,'Effective density for sh=',shapenoise,': '+ntostr(effective_density)

    ;
    ; Now the same but as a function of magnitude
    ; Mean weight as a function of magnitude
    ;

    magbin = 0.25
    bs = binner(struct.oldstruct[wgal].mag_auto, weights_scaled, $
        binsize=magbin, rev=rev)
    ;nperbin = 100
    ;bs = binner(struct.oldstruct[wgal].mag_auto, weights_scaled, $
    ;    nperbin=nperbin, rev=rev)

    ; plot weight vs mag
    xt = 'IAB'
    yt = 'weight'
    xrange = [struct.minmag,struct.maxmag]
    pplot, bs.xmean, bs.ymean, yerr=bs.yerr, psym=8, aspect=1, $
        yrange = [0.0, 1.2], xrange=xrange, /xsty, xtit=xt, ytit=yt
    pplot, $
        [bs.xmean, max(bs.xmean)+magbin/2.],  $
        [bs.hist/float(max(bs.hist))/2.,0.0],psym=10, $
        /overplot

    legend,'IAB < '+ntostr(struct.maxmag,f='(I0.2)'), box=0

    key=prompt_kbrd()

    nbin = n_elements(bs.xmean)
    mag_auto = bs.xmean
    mshear_sens = replicate(9999.0, nbin)
    meffective_density = replicate(0.0, nbin)
    for i=0L, nbin-1 do begin
        if rev[i] ne rev[i+1] then begin
            w=rev[ rev[i]:rev[i+1]-1 ]
            if n_elements(w) gt 1 then begin
                wmom, e1[w], blah, we1mean, we1sig, we1err, weights=weights[w]
                wmom, e2[w], blah, we2mean, we2sig, we2err, weights=weights[w]

                ; This is shear sensitivity for this field.
                mshear_sens[i] = 0.5*mean([we1err, we2err])
                meffective_density[i] = (shapenoise*0.5/mshear_sens[i])^2/area

            endif 
        endif
    endfor

    ; Adding new exposures
    ; I think this is wrong
    nexp3 = 3
  
    shear_sens3_squared = $
        shear_sens^2/nexp3 + shear_sens_shonly^2*(1.0 - 1.0/nexp3)
    shear_sens3 = sqrt( shear_sens3_squared )
    effdens3 = (shapenoise*0.5/shear_sens3)^2/area


    ; more correct way to estimate error with new exposures.  This calculates
    ; the expected relative error for n equivalent images compared to this
    ; one image

    mst = hgsc_calc_multi_factors(nexp, ellip_err, shapenoise)

    ; neffective stuff for all
    nstruct = {$
        seeing: struct.seeing, $
        med_seeing: struct.med_seeing, $
        shapenoise:shapenoise, $
        area:area, $
        $
        raw_density_unweighted: raw_density, $
        raw_density:weighted_density, $
        $
        shear_sens: shear_sens, $
        effective_density:effective_density, $
        $
        shear_sens_shonly: shear_sens_shonly, $
        effective_density_shonly: effective_density_shonly, $
        $
        shear_sens3: shear_sens3, $
        effective_density3: effdens3, $
        $
        multi_nexp: nexp, $
        multi_sens_fac: mst.multi_sens_fac, $
        multi_effdens_fac:mst.multi_effdens_fac, $
        $
        mag_auto:mag_auto, $
        mweight: bs.ymean ,$
        mweight_err: bs.yerr, $
        mhist: bs.hist, $
        mdensity:bs.hist/area, $
        mshear_sens:mshear_sens, $
        meffective_density:meffective_density $
    }

    return, nstruct
end

function huan_goods_sensitivity_convolved, seeing, $
        struct=struct, $
        wgood=wgood, $
        wstar=wstar, $
        wgal=wgal, $
        r=r, $
        type=type, $
        fwhm_range=fwhm_range,$
        use_cursor=use_cursor

    common _htsc_common, shapenoise, area, survey_area, nexp
    shapenoise=0.32
    area = 61.2958
    ; This is for dark energy survey
    survey_area = 5000.0
    ; Different equivalent exposures
    nexp = [1.0,1.5,2.0,2.5,3.0]


    ; for extracting galaxies
    minmag = 20.0
    maxmag = 30.0
    magstr = ntostr(maxmag,4)
    rcut = 0.9

    seestr = ntostr(seeing, 4, /round)+'_'
    if n_elements(type) eq 0 then typestr = '' else typestr=type+'_'
    dir = esheldon_config('des_sensitivity_dir')

    if n_elements(struct) eq 0 then begin 
        catdir = path_join(dir, 'cat')
        pattern = path_join(catdir, typestr+seestr+"*admom*.fit")
        files = findfile(pattern)
      
        struct = mrdfits_multi(files)
    endif 

    seeing_guess = sqrt(0.1^2 + seeing^2)
    newstruct = hgsc_get_gals($
        seeing_guess, struct, rcut, minmag, maxmag,$
        fwhm_range=fwhm_range,$
        use_cursor=use_cursor)

    key=prompt_kbrd()

    sens_struct = hgsc_calc_weights_and_dens(newstruct)
    return, sens_struct

end 


