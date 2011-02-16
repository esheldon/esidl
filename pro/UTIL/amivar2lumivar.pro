function amivar2lumivar, amivar, lum, err=err

    if n_elements(amivar) eq 0 OR n_elements(lum) eq 0 then begin 
        print,'-syntax: lumivar=amivar2lumivar(amivar, lum, err=)'
        print
        on_error, 2
        message,'Halting'
    endif 

    lum_ivar = amivar*0
    w=where(lum gt 0, nw)
    if nw ne 0 then begin 
        lum_ivar[w] = amivar[w]/lum[w]^2/(0.4*alog(10.))^2
    endif 

    if arg_present(err) then begin
        err = lum_ivar*0
        good=where(lum_ivar gt 0, ngood, comp=bad, ncomp=nbad)
        if ngood gt 0 then err[good] = sqrt(1.0/lum_ivar[good])
        if nbad gt 0 then err[bad] = 9999.0
    endif
    return, lum_ivar

end 



