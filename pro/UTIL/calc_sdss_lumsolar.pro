function _csl, absmag, clr

    sun=[6.38,5.06,4.64,4.53,4.52]

    ;; the 25 is to put in units of 10^10

    num=n_elements(absmag)
    lumsolar = replicate(-9999.0, num)

    w = where(absmag gt -30 and absmag lt -5, nw)
    if nw ne 0 then begin 
        lumsolar[w] = 10.^( (absmag[w] - sun[clr] + 25.0)/(-2.5) )
    endif 
    return,lumsolar

end

function _calc_sdss_absmag_from_lumsolar, lumsolar, clr
  
    sun=[6.38,5.06,4.64,4.53,4.52]

    num = n_elements(lumsolar)
    absmag = replicate(-9999.0, num)

    ;; The 25.0 assumes lumsolar in 10^10
    w=where(lumsolar GT 0 AND lumsolar LT 1.e4, nw)
    if nw ne 0 then begin 
        absmag[w] = sun[clr] - 25.0 - 2.5*alog10(lumsolar[w])
    endif 
    return,absmag

end 

function calc_sdss_lumsolar, absmag, $
        clr=clr, $
        amerr=amerr, lumerr=lumerr, $
        amivar=amivar, lumivar=lumivar, $
        inverse=inverse

    if n_params() lt 1 then begin 
        print,'-Syntax: lumsolar = calc_sdss_lumsolar(absmag, clr=, /inverse)'
        print
        print,'Returns luminosity in units of 10^10 solar'
        print,'If clr= not sent, assumes [5,N] array.  clr must be scalar'
        print,'/inverse: the argument is lumsolar (10^10 solar) and returned is a bsmag'
        print
        message,'Halting'
    endif 

    nClr = n_elements(clr)
    if (nclr ne 0) then begin 
        if (nclr ne 1) then message,'clr=clr should be a single bandpass'
    endif 

    if n_elements(amerr) ne 0 then begin
        amivar = amerr*0
        w=where(amerr gt 0, nw)
        if nw ne 0 then amivar[w] = 1.0/amerr^2
    endif

    nobj = n_elements(absmag[0, *])

    if nclr ne 0 then begin 
        if keyword_set(inverse) then begin
            return, _calc_sdss_absmag_from_lumsolar(absmag, clr)
        endif else begin 
            lumsolar=_csl(absmag, clr)
            if n_elements(amivar) ne 0 then begin
                lumivar=amivar2lumivar(amivar, lumsolar, err=lumerr)
            endif
            return, lumsolar
        endelse 
    endif else begin 

        ;; do all five bandpasses
        nclr = 5
        arr = replicate(-9999.0, 5, nobj)
        if n_elements(amivar) ne 0 then begin
            arr_ivar = arr*0
            if arg_present(
        endif
        for clr=0l, nclr-1 do begin 

            sendarr = reform( absmag[clr,*] )
            if keyword_set(inverse) then begin
                arr[clr,*] = _calc_sdss_absmag_from_lumsolar(sendarr, clr)
            endif else begin 
                arr[clr,*] = _csl(sendarr, clr)
            endelse 

        endfor 

        return,arr
    endelse 

end 



