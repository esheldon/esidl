function fitngauss, x, y, $
    error=error, $
    aguess=aguess, $
    parinfo=parinfo, $
    yfit=yfit, sigma=sigma, $
    silent=silent


    if n_params() lt 2 then begin 
        on_error, 2
        print,'-Syntax: parms = fitngauss(x, y, error=, aguess=, parinfo=, yfit=, sigma=, /silent)'
        print,'Either send aguess or parinfo'
        print,''
        message,'Halting'
    endif 

    npar = n_elements(parinfo)
    naguess = n_elements(aguess)
    if npar eq 0 and naguess eq 0 then begin
        on_error, 2
        message,'You must send parinfo or aguess'
    endif

    nx = n_elements(x)
    ny = n_elements(y)
    if ny ne nx then begin
        print,'X and Y must be same size'
        on_error, 2
        message,'Halting'
    endif 

    if n_elements(error) eq 0 then error = replicate(1d, nx)
    w=where(error EQ 0., nw)
    if nw ne 0 then begin
        print,'Errors must be non-zero'
        on_error, 2
        message,'Halting'
    endif 

    if npar eq 0 then begin

        if naguess mod 3 ne 0 then begin
            message,'npar must be a multiple of 3'
        endif
        ngauss = naguess/3

        parinfo = replicate({value:0.0, limited:[0,0], limits:[0d,0d]},naguess)
        parinfo.value = aguess

        for i=0L, ngauss-1 do begin

            ii = i*3
            ; amplitude must be reasonably large
            parinfo[ii].limited[0] = 1
            parinfo[ii].limits[0] = max(y)/100

            ; mean must be within the range
            parinfo[ii+1].limited = 1
            parinfo[ii+1].limits = [min(x), max(x)]

            ; sigma must be positive
            parinfo[ii+2].limited[0] = 1
            parinfo[ii+2].limits[0] = 1.e-2*(max(x)-min(x))
        endfor

    endif else begin
        if npar mod 3 ne 0 then begin
            messge,'npar must be a multiple of 3'
        endif
        ngauss = npar/3
    endelse

    itmax = 300 
    MYFUNCT = 'ngauss'
    a = MPFITFUN(MYFUNCT, x, y, error, parinfo=parinfo, $
                 AUTODERIVATIVE=1, MAXITER=itmax, niter=iter, $
                 yfit=yfit, perror=sigma, covar=covar,/quiet,/double,$
                 status=mpstatus)

    if not keyword_set(silent) then begin 

        chi2 = total( ( (y-yfit)/error )^2, /double)
        dof = n_elements(y)-( ngauss*3 )
        chi2per = chi2/dof
        print
        print,'--------------------------------------------------------'
        print,'best fit: '
        for i=0L, ngauss-1 do begin
            ii=i*3
            istr = ntostr(i+1)
            print,'   norm'+istr+': '+$
                ntostr(a[ii])+!plusminus+ntostr(sigma[ii])
            print,'   mean'+istr+': '+$
                ntostr(a[ii+1])+!plusminus+ntostr(sigma[ii+1])
            print,'   sig'+istr+':  '+$
                ntostr(a[ii+2])+!plusminus+ntostr(sigma[ii+2])
        endfor
        print,'iterations = ',iter
        print,'chi^2/dof  = ',ntostr(chi2)+'/'+ntostr(dof)+' = '+ntostr(chi2per)
    endif 

    return, a
end 
