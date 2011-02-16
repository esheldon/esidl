function lumsolar2am, lumsolar, sunmags

    if n_elements(lumsolar) eq 0 or n_elements(sunmags) eq 0 then begin
        print,'-Syntax: absmag=lumsolar2am(lumsolar, sunmags)'
        on_error, 2
        message,'Halting'
    endif
 
    absmag=-2.5*alog10(lumsolar) + sunmags
    return, absmag
end
