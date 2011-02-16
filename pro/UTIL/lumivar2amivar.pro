function lumivar2amivar, lum, lumivar, err=err

    if n_elements(lum) eq 0 or n_elements(lumivar) eq 0 then begin
        print,'-Syntax: amivar = lumivar2amivar(lum, lumivar, err=)'
        on_error, 2
        message,'Halting'
    endif

    amivar = lumivar*lum^2*(0.4*alog(10.))^2
    if arg_present(err) then begin
        err=amivar*0
        good = where(amivar gt 0, ngood, comp=bad, ncomp=nbad)
        if ngood gt 0 then err[good] = sqrt(1.0/amivar[good])
        if nbad gt 0 then err[bad] = 9999.0
    endif
    return, amivar
end
