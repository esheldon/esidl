function _fitm2l_mcmc_loglike, par
    common _m2l_block, r, m2l, cinv, psigma, cov_is_error

    ; simple priors
    ;if par[0] lt 0 or par[0] gt 5000 then return, -10000.0
    if par[0] lt 0 or par[0] gt 800 then return, -10000.0
    if par[1] lt 0 or par[1] gt 20 then return, -10000.0
    if par[2] lt 0 or par[2] gt 5 then return, -10000.0

    func = m2lfunc(r, par)
    diff = m2l-func

    if cov_is_error then begin
        chi2 = -total(cinv*diff^2)
    endif else begin
        chi2 = -transpose(diff)#( reform(cinv##diff) ) 
    endelse
    return, chi2
end

function _fitm2l_mcmc_step, seed, par
    common _m2l_block, r, m2l, cinv, psigma, cov_is_error

    npar = n_elements(psigma)
    return, par+psigma*randomn(seed, npar)
end

pro _fitm2l_mcmc_set_common, rin, m2lin, m2lcov, psigmain, cov_is_error_in

    common _m2l_block, r, m2l, cinv, psigma, cov_is_error

    cov_is_error = cov_is_error_in
    r = rin
    m2l = m2lin
    psigma = psigmain

    if cov_is_error then begin
        cinv = 1d/m2lcov^2
    endif else begin
        cinv = invert(m2lcov, status)
        
        if status ne 0 then begin
            message,'Failed to invert covariance matrix'
        endif
    endelse

 
end

function _fitm2l_mcmc_geterror, m2lcov, cov_is_error=cov_is_error
    ; Determine whether error or covariance matrix was sent
    dim_of_cov = size(m2lcov, /n_dim)
    if dim_of_cov eq 1 then begin
        error = m2lcov
        cov_is_error = 1
    endif else begin
        error = cov2diag(m2lcov)
        cov_is_error = 0
    endelse

    return, error
end
function fitm2l_mcmc, r, m2l, m2lcov, nstep, burnin, doplot=doplot

    ;; First fit without cov
    error = _fitm2l_mcmc_geterror(m2lcov, cov_is_error=cov_is_error)
    fitst = fitm2l(r, m2l, error)

    ; starting guess
    aguess = [fitst.m2l, fitst.rs, fitst.index]

    ; for stepping
    psigma = [$
        fitst.m2l_err > 0.05*abs(fitst.m2l) < abs(fitst.m2l), $
        fitst.rs_err > 0.05*abs(fitst.rs) < abs(fitst.rs), $
        fitst.index_err > 0.05*abs(fitst.index) < abs(fitst.index)]

    ; Now set common block
    _fitm2l_mcmc_set_common, r, m2l, m2lcov, psigma, cov_is_error

    mcmc=obj_new('mcmc')

    trials = mcmc->run('_fitm2l_mcmc_step', '_fitm2l_mcmc_loglike', $
                        nstep, aguess, $
                        like=like, /log)

    a = mcmc->extract_stats(trials, like, burnin)

    if keyword_set(doplot) then begin
        mcmc->plot1d, trials, like, burnin, names=[!csym.gamma_cap, 'r!Ds!N',!csym.alpha]
        mcmc->plot2d, trials, like, burnin, names=[!csym.gamma_cap, 'r!Ds!N',!csym.alpha]
    endif

    minchi2 = min(abs(like))
    dof = n_elements(r)-(3-1)
    chi2per = minchi2/dof

    yfit = m2lfunc(r, [ a[0,0], a[1,0], a[2,0] ] )
    fitst_mcmc = { $
        m2l: a[0,0],$
        m2l_err: a[0,1], $
        rs: a[1,0], $
        rs_err: a[1,1], $
        index: a[2,0], $
        index_err: a[2,1], $
        yfit:yfit, $
        chi2: minchi2, $
        dof: dof, $
        chi2reduced: chi2per, $
        $
        m2l_0: fitst.m2l, $
        m2l_err_0: fitst.m2l_err, $
        rs_0: fitst.rs, $
        rs_err_0: fitst.rs_err, $
        index_0: fitst.index, $
        index_err_0: fitst.index_err, $
        yfit_0: fitst.yfit, $
        chi2_0: fitst.chi2, $
        dof_0: fitst.dof, $
        chi2reduced_0: fitst.chi2reduced $
    }

    IF NOT keyword_set(silent) THEN BEGIN 

        print
        print,'--------------------------------------------------------'
        print,'best mcmc fit after burnin: '
        print,'  m2l: '+ntostr(a[0,0])+!plusminus+ntostr(a[0,1])
        print,'  rs: '+ntostr(a[1,0])+!plusminus+ntostr(a[1,1])
        print,'  index: '+ntostr(a[2,0])+!plusminus+ntostr(a[2,1])

        print,'  chi^2/dof  = ',ntostr(minchi2)+'/'+ntostr(dof)+' = '+ntostr(chi2per)
    ENDIF 

    return, fitst_mcmc 
end
