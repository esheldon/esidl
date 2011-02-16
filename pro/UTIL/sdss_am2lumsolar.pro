function sdss_am2lumsolar, absmag, $
        clr=clr, band_shift=band_shift, $
        amerr=amerr, lumerr=lumerr, $
        amivar=amivar, lumivar=lumivar

    if n_params() lt 1 then begin 
        print,'-Syntax: lumsolar = sdss_am2lumsolar(absmag, clr=, band_shift=, amerr=, lumerr=, amivar=, lumivar=)'
        print
        print,'Returns luminosity in units of solar'
        print,'If clr= not sent, assumes [5,N] array.  clr must be scalar'
        on_error, 2
        message,'Halting'
    endif 

    sunmags=sdss_solarmags(band_shift=band_shift)

    nClr = n_elements(clr)
    if (nclr ne 0) then begin 
        if (nclr ne 1) then message,'clr=clr should be a single bandpass'
    endif else begin
        sz=size(absmag)
        if sz[0] ne 0 or sz[1] ne 5 then begin
            message,'Unless clr= is sent, absmag must be a [5,N] array'
        endif
    endelse

    if n_elements(amerr) ne 0 then begin
        amivar = amerr*0
        w=where(amerr gt 0, nw)
        if nw ne 0 then amivar[w] = 1.0/amerr[w]^2
    endif

    nobj = n_elements(absmag[0, *])

    if nclr ne 0 then begin 
        lumsolar=am2lumsolar(absmag, sunmags[clr])
        if n_elements(amivar) ne 0 then begin
            lumivar=amivar2lumivar(amivar, lumsolar, err=lumerr)
        endif
    endif else begin 

        ;; do all five bandpasses
        nclr = 5
        lumsolar = replicate(-9999.0, 5, nobj)
        if n_elements(amivar) ne 0 then begin
            lumivar = lumsolar*0
            if arg_present(lumerr) then lerr=lumivar+9999
        endif
        for clr=0l, nclr-1 do begin 

            amtmp = reform( absmag[clr,*] )
            lumsolar[clr,*] = am2lumsolar(amtmp, sunmags[clr])

            if n_elements(amivar) ne 0 then begin
                lumivar[clr,*] = $
                    amivar2lumivar(amivar[clr,*], lumsolar[clr,*], err=err)

                if arg_present(lumerr) then lumerr[clr,*] = err
            endif

        endfor 
    endelse 

    return, lumsolar
end 



