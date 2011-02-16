;+
; NAME:
;	icl_fakeimage()
;
; PURPOSE:
;	create a fake image with icl, either a sersic profile or a double-dev
;	profile.  Uses parameters drawn from the Gonzalez et al. paper, which
;	are read from disk.  The location is hard coded at the top of the
;	program.  No noise is added.
;
; CALLING SEQUENCE:
;	im=icl_fakeimage(type, z, index=random, sblim=28, pixscale=0.396, 
;		omega_m=0.3, h=1.0, /write=, dir=
;
; INPUTS:
;	type: Either "sersic" or "2dev"
;	z: redshift
;
; OPTIONAL INPUTS:
;	index: The index in gonzales catalog.  This object will be used to
;		generate an ICL image.  Default is random.
;	sblim:  Surface brightness limit, default 28 mags/square arcsec.	
;	pixscale:  Default 0.396''/pixel, SDSS
;	omega_m: Default 0.3
;	h: default 1.0
;
; OUTPUTS:
;	An image with the fake ICL placed in.
;
; MODIFICATION HISTORY:
;	Documented: 2008-11-05, Erin Sheldon, BNL
;
;-

function icl_sim_dir
   return,'/global/early2/esheldon/icl/sim'
end
function icl_sim_file, type, num, dir=dir
    dir=icl_sim_dir()
    file = 'icl-sim-'+strlowcase(type)+ntostr(num,format='(I03)')+'.fits'
    file=concat_dir(dir,file)
    return, file
end



function icl_2dev, r50_1, r50_2, siz, $
        aratio1=aratio1, theta1=theta1, counts1=counts1, core1=core1, $
        aratio2=aratio2, theta2=theta2, counts2=counts2, core2=core2, $
        blanton=blanton

    if n_elements(r50_1) eq 0 and n_elements(r50_2) eq 0 then begin
        print,'-Sytnax: icl_2dev(r50_1, r50_2, [xsize,ysize], '
        print,'                  counts1=, aratio1=, theta1=, core1=, '
        print,'                  counts2=, aratio2=, theta2=, core2='
        on_error, 2
        message,'Halting'
    endif

    if keyword_set(blanton) then begin
        im1 = dfakegal(flux=counts1, sersicn=4.0, r50=r50_1, nx=siz[0], ny=siz[1]) 
        im2 = dfakegal(flux=counts2, sersicn=4.0, r50=r50_2, nx=siz[0], ny=siz[1]) 
    endif else begin
        im1 = sersic_image(r50_1, 4.0, siz, counts=counts1, aratio=aratio1, theta=theta1, core=core1)
        im2 = sersic_image(r50_2, 4.0, siz, counts=counts2, aratio=aratio2, theta=theta2, core=core2)
    endelse

    im = im1+im2
    return,im
end


function icl_distmod,z, h=h, omega_m=omega_m
    dmpc = angdist_lambda(z, h=h, omegamat=omega_m)
    dpc = dmpc*1.e6
    distmod = 5.0*alog10(dpc/10.0)
    return, distmod
end

function icl_absmag2nmgy, absmag, z, h=h, omega_m=omega_m, mag=mag
    if n_elements(absmag) eq 0 or n_elements(z) eq 0 then begin
        print,'-Syntax: nmgy = icl_absmag2nmgy(absmag, z, h=, omega_m=, mag=)'
        on_error, 2
        message,'Halting'
    endif
    ; convert absmag to apparrent mag
    distmod = icl_distmod(z, h=h, omega_m=omega_m)

    mag=absmag + distmod

    ; mag = 22.5 - 2.5*alog10(nmgy)
    nmgy = 10.0^( -0.4*(mag-22.5) )
    return, nmgy
end

; convert a magnitudes per square arcsec type of number into
; the actual flux density in nanomaggies per square arcsec.
;  sb (mag/arcsec^2) = 22.5 - 2.5*alog10( flux/arcsec^2 )
;  thus flux density = 10^( -0.4*(sb-22.5) )
function icl_sb2fluxdensity, sb
    fd = 10.0^(-0.4*(sb-22.5) ) 
    return, fd
end

function icl_2sersic_sbrad, A1, r50_1, n1, A2, r50_2, n2, sblimit
    ; this will give half the flux density
    sbuse = sblimit + 2.5*alog10(2.0)
    r1 = sersic_sbrad(A1, r50_1, n1, sbuse)
    r2 = sersic_sbrad(A2, r50_2, n2, sbuse)
    return, max([r1,r2])
end


function icl_sizeconvert, z, rkpc, omega_m=omega_m, h=h, pixscale=pixscale
    if n_elements(z) eq 0 or n_elements(rkpc) eq 0 then begin
        print,'-Syntax: st=icl_sizeconvert(z, rkpc, pixscale=)'
        on_error, 2
        message,'Halting'
    endif

    DA = angdist_lambda(z, h=h, omegamat=omega_m)

    rmpc = rkpc/1000.0
    rradians = rmpc/DA
    rdegrees = rradians*180.0/!pi
    rarcsec = rdegrees*3600.0

    st = {kpc:rkpc, mpc:rmpc, radians:rradians, degrees:rdegrees, arcsec:rarcsec}
    if arg_present(pixscale) then begin
        rpixels = rarcsec/pixscale
        st=create_struct(st, 'pixels', rpixels)
    endif
    return, st
end

function icl_fakeimage, type, z, index=index, sblim=sblim, pixscale=pixscale, omega_m=omega_m, h=h, core=core, write=write, dir=dir

    if n_elements(type) eq 0 or n_elements(z) eq 0 then begin
        print,'-Syntax: image = icl_fakeimage(type, z, index=random, sblim=28, pixscale=0.396, omega_m=0.3, h=1.0, /write, dir=)'
        on_error, 2
        message,'Halting'
    endif

    ; arcsec per pixel
    if n_elements(pixscale) eq 0 then pixscale = 0.396 
    if n_elements(sblim) eq 0 then sblim=28.0
    if n_elements(omega_m) eq 0 then omega_m = 0.3
    if n_elements(h) eq 0 then h=1.0
    if n_elements(core) eq 0 then core = 0.3 ; kpc. worked for index=3, A2400

    ; read gonzales data 
    struct = icl_gonzales_read(type)
    nst = n_elements(struct)

    ; get a random index if not sent
    if n_elements(index) eq 0 then begin
        rr=randomu(seed, nst)
        s = sort(rr)
        index = s[0]
    endif

    object = struct[index]

    case strlowcase(type) of
        'sersic': begin
            ; get radii in various units
            r50st = icl_sizeconvert(z, object.r50, omega_m=omega_m, h=h, $
				pixscale=pixscale)
            cst = icl_sizeconvert(z, core, omega_m=omega_m, h=h, $
				pixscale=pixscale)

            ; convert absolute mag to a flux at given redshift
            nmgy = icl_absmag2nmgy(absmag, z, h=h, omega_m=omega_m, mag=mag)

            ; convert surface brightness limit to a flux density limit
            fdlim = icl_sb2fluxdensity(sblim)

            ; get radius at which profile falls to input surface brightness
            maxrad_arcsec = sersic_sbrad(nmgy, r50st.arcsec, object.n, sblim)
            npix = long(2*maxrad_arcsec/pixscale)

            image = sersic_image(r50st.pixels, object.n, $
				[npix,npix], counts=nmgy, core=cst.pixels)
            return, image
        end
        '2dev': begin
            r50_1st = icl_sizeconvert(z, struct[index].r50_1, $
				omega_m=omega_m, h=h, pixscale=pixscale) 
            r50_2st = icl_sizeconvert(z, struct[index].r50_2, $
				omega_m=omega_m, h=h, pixscale=pixscale) 

            ; convert absolute mag to a flux at given redshift
            nmgy1 = icl_absmag2nmgy(struct[index].absmag1, z, h=h, $
				omega_m=omega_m, mag=mag1)
            nmgy2 = icl_absmag2nmgy(struct[index].absmag2, z, h=h, $
				omega_m=omega_m, mag=mag2)

            ; convert surface brightness limit to a flux density limit
            fdlim = icl_sb2fluxdensity(sblim)

            ; get radius at which profile falls to input surface brightness
            maxrad_arcsec = icl_2sersic_sbrad(nmgy1, r50_1st.arcsec, 4.0, $
				nmgy2, r50_2st.arcsec, 4.0, sblim)
            npix = long(2*maxrad_arcsec/pixscale)

            dev2 = icl_2dev(r50_1st.pixels, r50_2st.pixels, [npix,npix], $
				counts1=nmgy1, counts2=nmgy2, /blanton)
            return, dev2
        end
        else: message,'Bad type: '+ntostr(type)
    endcase

    if keyword_set(write) then begin
        outfile = icl_sim_file(type,index, dir=dir)
        print
        print,'Writing to file: ',outfile
        writefits, outfile, image

        hdr=headfits(outfile)
        case strlowcase(type) of
            'sersic': begin
                sxaddpar, hdr, 'name', ntostr(object.name)
                sxaddpar, hdr, 'z', z
                sxaddpar, hdr, 'r50kpc', r50st.kpc
                sxaddpar, hdr, 'core_kpc', cst.kpc
                sxaddpar, hdr, 'absmag', absmag
                sxaddpar, hdr, 'mag', mag
                sxaddpar, hdr, 'nmgy', nmgy
                sxaddpar, hdr, 'sblim', sblim
                sxaddpar, hdr, 'fdlim', fdlim
                sxaddpar, hdr, 'nx', npix
                sxaddpar, hdr, 'ny', npix
                modfits, outfile, 0, hdr
            end
            '2dev': begin
                sxaddpar, hdr, 'name', ntostr(object.name)
                sxaddpar, hdr, 'z', z
                sxaddpar, hdr, 'r50kpc_1', r50st_1.kpc
                sxaddpar, hdr, 'r50kpc_2', r50st_2.kpc
                sxaddpar, hdr, 'core_kpc', cst.kpc
                sxaddpar, hdr, 'absmag1', struct[index].absmag1
                sxaddpar, hdr, 'absmag2', struct[index].absmag2
                sxaddpar, hdr, 'mag1', mag1
                sxaddpar, hdr, 'mag2', mag2
                sxaddpar, hdr, 'nmgy1', nmgy1
                sxaddpar, hdr, 'nmgy2', nmgy2
                sxaddpar, hdr, 'sblim', sblim
                sxaddpar, hdr, 'fdlim', fdlim
                sxaddpar, hdr, 'nx', npix
                sxaddpar, hdr, 'ny', npix
                modfits, outfile, 0, hdr
            end
            else: message,'Bad type: '+ntostr(type)
        endcase
    endif
    !p.multi=0

    return, image
end

