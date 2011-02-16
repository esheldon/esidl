function sdss_lumsolar2am, lumsolar, $
        clr=clr, band_shift=band_shift, $
        amerr=amerr, lumerr=lumerr, $
        amivar=amivar, lumivar=lumivar

    if n_params() lt 1 then begin 
        print,'-Syntax: absmag = sdss_lumsolar2absmag(lumsolar, clr=, band_shift=, lumerr=, amerr=, lumivar=, amivar=)'
        print
        print,'lumsolar is luminosity in units of solar. Returns absolute'
        print,'  magnitude'
        print,'If clr= not sent, assumes [5,N] array.  clr must be scalar'
        on_error, 2
        message,'Halting'
    endif 

    sunmags=sdss_solarmags(band_shift=band_shift)

    nClr = n_elements(clr)
    if (nclr ne 0) then begin 
        if (nclr ne 1) then message,'clr=clr should be a single bandpass'
    endif else begin
        sz=size(lumsolar)
        if sz[0] ne 0 or sz[1] ne 5 then begin
            message,'Unless clr= is sent, lumsolar must be a [5,N] array'
        endif
    endelse

    if n_elements(lumerr) ne 0 then begin
        lumivar = lumerr*0
        w=where(lumerr gt 0, nw)
        if nw ne 0 then lumivar[w] = 1.0/lumerr[w]^2
    endif

    nobj = n_elements(lumsolar[0, *])

    if nclr ne 0 then begin 
        absmag=lumsolar2am(lumsolar, sunmags[clr])
        if n_elements(lumivar) ne 0 then begin
            amivar=lumivar2amivar(lumsolar, lumivar, err=amerr)
        endif
    endif else begin 

        ;; do all five bandpasses
        nclr = 5
        absmag = replicate(-9999.0, 5, nobj)
        if n_elements(lumivar) ne 0 then begin
            amivar = lumsolar*0
            if arg_present(amerr) then amerr=amivar+9999
        endif
        for clr=0l, nclr-1 do begin 

            ltmp = reform( lumsolar[clr,*] )
            absmag[clr,*] = lumsolar2am(ltmp, sunmags[clr])

            if n_elements(lumivar) ne 0 then begin
                amivar[clr,*] = $
                    lumivar2amivar(lumsolar[clr,*], lumivar[clr,*], err=err)

                if arg_present(amerr) then amerr[clr,*] = err
            endif

        endfor 
    endelse 

    return, absmag
end 



