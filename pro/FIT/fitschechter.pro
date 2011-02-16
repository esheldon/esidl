function fitschechter, x, y, yerror, $
    aguess=aguess, $
    params=params, $
    sigma=sigma, $
    parinfo=parinfo, $
    linear=linear, $
    yfit=yfit, $
    silent=silent


    if n_params() lt 2 then begin 
        on_error, 2
        print,'-Syntax: fitstruct = fitschechter(x, y, yerror, aguess=, parinfo=, /linear, params=, sigma=, yfit=, /silent)'
        print,'If aguess/parinfo not sent, a guess is made'
        print,''
        message,'Halting'
    endif 

    npar = n_elements(parinfo)
    naguess = n_elements(aguess)

    if npar eq 0 and naguess eq 0 then begin
        aguess = [1.d-2, -20d, -1d]
    endif
    ; very basic limits
    if npar eq 0 then begin
        parinfo = replicate({value:0.0, limited:[0,0], limits:[0d,0d]},3) 

        parinfo.value = aguess
        parinfo[0].limited[0] = 1
        parinfo[0].limits[0] = max(y)/100

        parinfo[2].limited = 1
        parinfo[2].limits = [-2d,2d]
    endif

    nx = n_elements(x)
    ny = n_elements(y)
    if ny ne nx then begin
        print,'X and Y must be same size'
        on_error, 2
        message,'Halting'
    endif 

    w=where(yerror EQ 0., nw)
    if nw ne 0 then begin
        print,'Errors must be non-zero'
        on_error, 2
        message,'Halting'
    endif 

    itmax = 300 
    if keyword_set(linear) then begin
        MYFUNCT = 'mpfit_schechter'
    endif else begin
        myfunct = 'mpfit_schechter_absmag'
    endelse

    a = MPFITFUN(MYFUNCT, x, y, yerror, parinfo=parinfo, $
                 AUTODERIVATIVE=1, MAXITER=itmax, niter=iter, $
                 yfit=yfit, perror=sigma, covar=covar,/quiet,/double,$
                 status=mpstatus)

    if not keyword_set(silent) then begin 

        chi2 = total( ( (y-yfit)/yerror )^2, /double)
        dof = n_elements(y)-3
        chi2per = chi2/dof
        print
        print,'--------------------------------------------------------'
        print,'best fit: '
        print,'  phistar: '+ntostr(a[0])+!plusminus+ntostr(sigma[0])
        print,'  mstar:   '+ntostr(a[1])+!plusminus+ntostr(sigma[1])
        print,'  alpha:   '+ntostr(a[2])+!plusminus+ntostr(sigma[2])
        print,'iterations = ',iter
        print,'chi^2/dof  = ',ntostr(chi2)+'/'+ntostr(dof)+' = '+ntostr(chi2per)
    endif 

    params=a
    fitstruct = $
        {phistar: a[0], phistar_err: sigma[0], $
         mstar: a[1], mstar_err: sigma[1], $
         alpha: a[2], alpha_err: sigma[2], $
         covar: covar, $
         chi2: chi2, $
         dof: dof, $
         chi2per: chi2per}
    return, fitstruct
end 
