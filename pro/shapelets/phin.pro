FUNCTION phin, n1, n2, x1, x2, scale

; This function returns the basis functions phi_n(x) for the hermite
; polynomials.  It is linked to code written by A. Refregier.
; 
; Author:  Brandon Kelly May 2002
;
; Input:
; n1 = first order indicator
; n2 = second order indicator
; x1 = array containing the values for x_1
; x2 = array containing the values for x_2
; scale= a two-element vector containing the x and y scales
;
; Output: the basis functions phi_n(x)

s1 = size(x1)
s2 = size(x2)
phi_n1 = dblarr(s1[1])
h_n1 = hermite(n1, x1 / scale[0])
phi_n1 = double( (2.0^n1 * sqrt(!pi) * factorial(n1) * scale[0])^(-0.5) * h_n1 * $
  exp( -(x1 / scale[0])^2 / 2.0) )

phi_n2 = dblarr(s2[1])
h_n2 = hermite(n2, x2 / scale[1])
phi_n2 = double((2.0^n2 * sqrt(!pi) * factorial(n2) * scale[1])^(-0.5) * h_n2 * $
  exp( -(x2 / scale[1])^2 / 2.0) )

phi = transpose(phi_n2) ## phi_n1

return, phi
end
