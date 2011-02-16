function errstring_expstring, expon, screen=screen
    if keyword_set(screen) then begin
        return, ntostr(expon)
    endif else begin
        return, '{'+ntostr(expon)+'}'
    endelse
end
function errstring_find_exponent, x
    return, long( alog10( abs(x) ) )
end

function errstring, val, error, valstring, errstring, screen=screen
    ; re-normalize such that the error is
    ; between 1 and 10

    err_expon = errstring_find_exponent(error)
    val_expon = errstring_find_exponent(val)

    tval = val/10d^err_expon
    terr = error/10d^err_expon

    tval = rnd(tval, 1)
    terr = rnd(terr, 1)

    if keyword_set(screen) then begin
        pm = string(177b)
        times = 'X'
    endif else begin
        pm = '\pm'
        times = '\times'
    endelse

    if err_expon eq -1 then begin
        f = '(F0.2)'

        tval = tval*10d^err_expon
        terr = terr*10d^err_expon
        valstring = ntostr(tval, f=f)
        errstring = ntostr(terr, f=f)
        pstring = valstring+pm+errstring
    endif else if err_expon eq 0 or err_expon eq 1 then begin
        f='(F0.1)'  
        tval = tval*10d^err_expon
        terr = terr*10d^err_expon
        valstring = ntostr(tval, f=f)
        errstring = ntostr(terr, f=f)
        pstring = valstring+pm+errstring
    endif else begin
        f='(F0.1)'  
        print,err_expon
        if err_expon gt 1 then begin
            tval=tval/10d
            terr=terr/10d
            err_expon = err_expon+1
        endif
        valstring = ntostr(tval, f=f)
        errstring = ntostr(terr, f=f)
        pstring = '('+valstring+pm+errstring+')'
        expstr = errstring_expstring(err_expon, screen=screen)
        pstring = pstring + times+'10^'+expstr
    endelse

    if not keyword_set(screen) then begin
        pstring = textoidl(pstring)
    endif
    return, pstring
end
