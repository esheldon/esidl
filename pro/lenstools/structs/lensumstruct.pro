function lensumstruct, arrval

    if n_params() lt 1 then begin 
        print,'-Syntax: lensum = lensumstruct(arrval)'
        return, -1
    endif 

    tmparrval = arrval
    tmparrval[*] = 1.e7

    lensumstruct = create_struct('zindex', 0L, $
                               'index', 0L, $
                               'totpairs', 0L, $
                               'ie', 1.0, $
                               'scritinv', 0., $
                               'DL', 0.0, $
                               'angMax', 0.0, $
                               'weight', 0.0, $
                               'pixelmaskflags', 0, $
                               $
                               'angsum', arrval, $
                               'rsum', arrval, $
                               'rmax_act', arrval, $
                               'rmin_act', tmparrval,$
                               $
                               'sigma', arrval, $
                               'sigmaerr', arrval, $
                               'sigerrsum', arrval, $ ;alternative error
                               'orthosig', arrval, $
                               'orthosigerr', arrval, $
                               'orthosigerrsum', arrval, $
                               $
                               'sshsum', 0., $
                               'wsum', arrval, $
                               'owsum', arrval, $
                               'wsum_ssh', 0., $
                               'npair', long(arrval))

    return, lensumstruct

end 

