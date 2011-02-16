function fitm2l, x, y, error, doplot=doplot

    ; asymptotic m/l, rs, slope
    aguess = [350.0d, 0.3d, 1.0d]

    func = 'm2lfunc'

    a = mpfitfun(func, x, y, error, aguess, $
        autoderivative=1, maxiter=300, niter=niter,$
        yfit=yfit, perror=sigma, covar=covar, /quiet, /double)

    IF NOT keyword_set(silent) THEN BEGIN 

        chi2 = total( ( (y-yfit)/error )^2, /double)
        dof = n_elements(y)-(3-1)
        chi2per = chi2/dof
        print
        print,'--------------------------------------------------------'
        print,'best fit: '
        print,'  m2l: '+ntostr(a[0])+!plusminus+ntostr(sigma[0])
        print,'  rs: '+ntostr(a[1])+!plusminus+ntostr(sigma[1])
        print,'  index: '+ntostr(a[2])+!plusminus+ntostr(sigma[2])

        print,'iterations = ',niter
        print,'chi^2/dof  = ',ntostr(chi2)+'/'+ntostr(dof)+' = '+ntostr(chi2per)
    ENDIF 


    if keyword_set(doplot) then begin
        yrange=prange(y, error, /slack)
        if yrange[0] lt 0 then begin
            yrange[0] = 1
        endif
        pplot, x, y, yerr=error, $
            psym=8, /xlog, /ylog, yrange=yrange, ystyle=3, aspect=1
        add_labels, ytickv=[200,300,400]
        add_labels, yaxis=1, ytickv=[200,300,400]
        oplot, x, yfit
    endif

    struct = {$
        x: x, y:y, error:error, $
        yfit: yfit, $
        chi2: chi2, $
        dof: dof, $
        chi2reduced: chi2per, $
        m2l: a[0], $
        m2l_err: sigma[0], $
        rs: a[1], $
        rs_err: sigma[1], $
        index: a[2], $
        index_err: sigma[2]}


    return, struct

end
