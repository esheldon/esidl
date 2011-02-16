pro linear_xi_sample

    log_rmin = -2.0
    log_rmax = 2.3
    r = 10d^( arrscl(dindgen(1000), log_rmin, log_rmax) )
    nplk = 200.0
    linear_xi, r, xi, $
        ns=0.98,$
        omega_m=0.25,$
        omega_b=0.055,$
        h=1.0,$
        sigma_8=1.0, $
        k=k, $
        nonorm=nonorm_pk

    f='~/tmp/pkidl.dat'
    print,'Writing to: ',f
    colprint, k, nonorm_pk, file=f

    f='~/tmp/xiidl.dat'
    print,'Writing to: ',f
    colprint, r, xi, file=f
end
