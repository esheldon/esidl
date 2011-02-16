function corr2cov, corr, diagerr

    ncorr=n_elements(corr) & nd=n_elements(diagerr)
    if ncorr eq 0 or nd eq 0 then begin 
        print,'-Syntax: cov = corr2cov(corr, diagerr)'
        on_error, 2
        message,'Halting'
    endif 

    ;; takes in a correlation matrix and diagonal errors and outputs the
    ;; corresponding covariance matrix
    ;; cov[i,u] = corr[i,j]*diag[i]*diag[j]

    siz = size(corr)

    if siz[0] ne 2 then message,'corr must be a matrix'
    nx = siz[1]
    ny = siz[2]

    if nx ne ny then message,'corr must be square'
    if nx ne nd then message,'diagerr must be same dimension in 1d as corr'

    cov = replicate(corr[0,0], nx, ny)
    cov[*] = 0.0

    for ix=0l, nx-1 do begin 
        for iy=0l, ny-1 do begin 
            cov[ix,iy] = corr[ix,iy]*diagerr[ix]*diagerr[iy]
        endfor 
    endfor 
  
    return, cov

end 
