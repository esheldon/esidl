pro _sigmacritinv_syntax
    on_error, 2 
    message,$
     '-Syntax: sc = scinv(zl, zs, omega_m=, omega_l=, omega_k=, flat=, h=, /cgs)', /inf
    message,'Possible inputs: ', /inf
    message,'  zl scalar, zs scalar', /inf
    message,'  zl scalar, zs vector', /inf
    message,'  zl vector, zs scalar', /inf
    message,'  zl vector, zs vector of the same length'
end


function sigmacritinv, zl, zs, $
        omega_m=omega_m, omega_l=omega_l, omega_k=omega_k, flat=flat, h=h, $
        cgs=cgs


    if n_params() lt 1 then begin
        print,'Syntax:  s = sigmacritinv(zlens, zsource, omega_m=, omega_l=, omega_k=, flat=, h=, /cgs)'
        print,'Returns the inverse lensing critical density in pc^2/msolar unless /cgs then cm^2/gm'
        on_error, 2
        message,'Halting'
    endif 

    nl = n_elements(zl) & ns = n_elements(zs)
    if nl eq 0 or ns eq 0 then begin
        _sigmacritinv_syntax
    endif
    ; bad case is where nl is not 1 and ns is not 1 but they are not equal
    if (nl ne 1) and (ns ne 1) then begin
        if nl ne ns then begin
            _sigmacritinv_syntax
        endif
    endif

    dl = angdist(0.0, zl, $
        omega_m=omega_m, omega_l=omega_l, omega_k=omega_k, flat=flat, h=h)
    ds = angdist(0.0, zs, $
        omega_m=omega_m, omega_l=omega_l, omega_k=omega_k, flat=flat, h=h)
    dls = angdist(zl, zs, $
        omega_m=omega_m, omega_l=omega_l, omega_k=omega_k, flat=flat, h=h)

    D = dls*dl/ds   ; Mpc/h
    D = D/1000.     ; Gpc/h

    if not keyword_set(cgs) then begin
        scinv = D/1663.0    ; (msolar/pc^2)^(-1)
    endif else begin
        scinv = D/0.3474    ; (gm/cm^2)^(-1)
    endelse

    ; for zsource < zlens we must return zero
    w=where(zs le zl, nw)
    if nw ne 0 then scinv[w] = 0.0

    return, scinv
end
