function am2lumsolar, absmag, sunmags

    if n_elements(absmag) eq 0 or n_elements(sunmags) eq 0 then begin
        print,'-Syntax: lumsolar=am2lumsolar(absmag, sunmags)'
        on_error, 2
        message,'Halting'
    endif
    arg = (absmag - sunmags)/(-2.5)
    lumsolar = 10d^arg
    return, lumsolar

end


