function chilevels, nu, perc=perc

    if n_elements(nu) eq 0 then begin
        print,'Usage:  lev = chilevels(nu, perc=)'
        print,' nu = # of parameters'
        on_error, 2
        message,'Halting'
    endif

    ; levels corresponding to 
    perc = [68.3, 90.0, 95.4, 99.0, 99.73, 99.99] ; % confidence levels
    case nu of
        1: lev = [1.00, 2.71, 4.00, 6.63, 9.00, 15.1]
        2: lev = [2.30, 4.61, 6.17, 9.21, 11.8, 18.4]
        3: lev = [3.53, 6.25, 8.02, 11.3, 14.2, 21.1]
        4: lev = [4.72, 7.78, 9.70, 13.3, 16.3, 23.5]
        5: lev = [5.89, 9.24, 11.3, 15.1, 18.2, 25.7]
        6: lev = [7.04, 10.6, 12.8, 16.8, 20.1, 27.8]
        else: message,'nu = 1 to 6 are supported'
    endcase

    return, lev
end
