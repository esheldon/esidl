;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Function to compute the inter-quartile range of a data set. This can
; give a more robust estimation of the standard deviation of a data
; set that is assumed to be approximately normal, as
;              
;                    IQR / 1.34 = Sigma
;                 
; Author : Brandon C. Kelly, Steward Obs., Sept. 2004
;
; INPUTS :
;    X - The data.
;
; OUTPUT :
;    The inter-quartile range of the data, IQR.
;
; OPTIONAL OUTPUT :
;    SIGMA - This will be a robust estimate of the standard deviation
;            assuming that the data set is drawn from a normal
;            distribution. SIGMA will be the smaller of the sample
;            standard deviation and IQR / 1.34.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function iqr, x, sigma

if n_params() eq 0 then begin
    print, 'Syntax- Result = iqr( x, [sigma])'
    return, 0
endif

nx = n_elements(x)

sorted = sort(x)

iqr = x[sorted[3 * nx / 4]] - x[sorted[nx / 4]]

sigma = min( [stddev(x, /nan), iqr / 1.34] )

return, iqr
end
