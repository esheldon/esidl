
function var, x, double = double, nan = nan, error=error

  on_error, 2

  error = 0.0
  nx = n_elements(x)
  if nx lt 2 then return,0. else begin
      result = moment( x, double=double, maxmoment=2, nan = nan)
      xvar = result[1]


      ; Calculate standard error on the variance
      if keyword_set(nan) then begin
          w = where( finite(x), nx)
          if nx eq 0 then begin
              return, 0
          endif
      endif
      error = sqrt(2.0/(nx-1))*xvar



      return, xvar
  endelse 
end
