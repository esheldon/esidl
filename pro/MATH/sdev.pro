
function sdev, x, double = double, nan = nan, error=error

  on_error, 2

  error = 0.0
  nx = n_elements(x)
  if nx lt 2 then return,0. else begin
      result = moment( x, double=double, maxmoment=2, nan = nan, sdev=tsdev)


      ; Calculate standard error on the standard deviation
      if keyword_set(nan) then begin
          w = where( finite(x), nx)
          if nx eq 0 then begin
              return, 0
          endif
      endif
      error = tsdev/sqrt(2*nx)



      return, tsdev
  endelse 
end
