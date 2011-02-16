function hermitecof, n
  
; July 99 - Written by A. Refregier
;
; PURPOSE: compute the polynomial coefficients for Hn(x), the Hermite
; polynomial of order n
; INPUT: n: order of the Hermite polynomial
; OUTPUT: hermitecof: coefficients ci for Hn(x)=Sum_i ci*x^i, i=0,..,n
  
c = dblarr(n+1)
for s=0, n/2 do $ 
  c(n-2*s) = (-1)^s*2.^(n-2*s)*factorial(n)/(factorial(s)*factorial(n-2*s))

return, c
end
  

