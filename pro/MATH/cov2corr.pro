function cov2corr, cov, fix=fix, status=status

    ncov = n_elements(cov)
    if ncov eq 0 then begin 
        on_error, 2
        print,'-Syntax: corr = cov2corr(cov, /fix, status=)'
        print,' status = 0 is good, 1 is negative Cii*Cjj, 2 is zero Cii*Cjj'
        message,'Halting'
    endif 

    status = 1

    ;; takes in a covariance matrix and outputs the
    ;; corresponding correlation matrix
    ;; corr[i,j] = cov[i,j]/sqrt(cov[i,i]*cov[j,j])

    siz = size(cov)

    if siz[0] ne 2 then message,'corr must be a matrix'
    nx = siz[1]
    ny = siz[2]

    if nx ne ny then message,'corr must be square'

    corr = replicate(cov[0,0], nx, ny)
    corr[*] = 0.0


    for ix=0l, nx-1 do begin 
        for iy=0l, ny-1 do begin 

            normfac = cov[ix,ix]*cov[iy,iy]

            if normfac gt 0.0 then begin 
                corr[ix,iy] = cov[ix,iy]/sqrt(normfac)
            endif else begin 

                if normfac lt 0.0 then begin                   
                    if keyword_set(fix) then begin 
                        corr[ix,iy] = 100.0
                    endif else begin 
                        message,'normalization term is negative',/inf
                        status = 1
                        return, -1
                    endelse 
                endif 
                if normfac eq 0.0 then begin

                    if keyword_set(fix) then begin 
                        corr[ix,iy] = 100.0
                    endif else begin 
                        message,'normalization term is zero',/inf
                        status = 2
                        return, -1
                    endelse 
                endif 

            endelse 
          
        endfor 
    endfor 
  
    status = 0
    return, corr

end 
