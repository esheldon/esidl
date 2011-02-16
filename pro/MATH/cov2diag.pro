function cov2diag, cov

    IF n_params() LT 1 THEN BEGIN 
        print,'-Syntax: diagerr = cov2diag(covariance)'
        return, -1
    ENDIF 

    sz = size(cov)

    nbin = sz[1]

    diagerr = dblarr(nbin)

    FOR i=0L, nbin-1 DO diagerr[i] = sqrt(cov[i, i])

    return, diagerr

end


