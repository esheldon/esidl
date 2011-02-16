function ngauss, x, p

    numgauss = n_elements(p)/3 

    model = double(x)
    model[*] = 0.0

    for i=0L, numgauss-1 do begin
        ii = i*3

        model[*] = model[*] + p[ii]*exp( -(x-p[ii+1])^2/2d/p[ii+2]^2 )

    endfor

    return, model
end 
