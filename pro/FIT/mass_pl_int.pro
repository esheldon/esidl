function mass_int_kernal1,x
return,(1d +x^2)^(!mass_int_kernal_exp) * (1d +2*x^2)
end

function mass_int_kernal2,x
;use taylor expansion for large numbers
if x lt 10d then  return,(1d +x^2)^(!mass_int_kernal_exp) * (1d +2*x*(x-sqrt(1d +x^2)))
if x ge 10d and x lt 100d then return,(1d +x^2)^(!mass_int_kernal_exp) * (1d/(4d*x^2)-1d/(8d*x^4)+5d/(16d*x^6))
if x gt 100d then return,(1d +x^2)^(!mass_int_kernal_exp) * (1d/(4d*x^2)-1d/(8d*x^4))
end

function mass_pl_int,X0,X1,b,out=out,in=in
;perform the numerical integral 
; \int_X0^X1 dx (1+x^2)^(-b/2) f(x)
;where f(x)=1+2*x^2 if "in" is set
; or f(x)=1+2*x^2 -2x*sqrt(1+x^2) in "out" is set

defsysv,"!mass_int_kernal_exp",-b/2.0
if keyword_set(in) then func="mass_int_kernal1" 
if keyword_set(out) then func="mass_int_kernal2" 
return,qromb(func,x0,x1,/double,eps=1e-8)
end


