function sigmacrit, zl, zs, $
        omega_m=omega_m, omega_l=omega_l, omega_k=omega_k, flat=flat, h=h, $
        cgs=cgs


    if n_params() lt 1 then begin
        print,'Syntax:  s = sigmacrit(zlens, zsource, omega_m=, omega_l=, omega_k=, flat=, h=, /cgs)'
        print,'Returns lensing critical density in msolar/pc^2 unless /cgs then gm/cm^2'
        on_error, 2
        message,'Halting'
    endif 

    scinv = sigmacritinv(zl, zs, $
        omega_m=omega_m, omega_l=omega_l, omega_k=omega_k, flat=flat, h=h, $
        cgs=cgs)
    return, 1/scinv


end 
