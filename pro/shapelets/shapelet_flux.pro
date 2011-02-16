;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Function to compute the flux of a galaxy from its shapelet
; coefficients.
;
; Author : B.C. Kelly, Steward Obs., Feb 2005
;
; Inputs :
;   COEFS - The shapelet coefficients.
;   SCALE - The scale of the shapelet functions.
;   N1, N2 - The N1 and N2 arrays that contain the order of the
;            shapelet coefficients. More specifically,
;            (N1[j],N2[j]) is the order of the jth shapelet. These
;            vectors should be the same size as COEFS.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function shapelet_flux, coefs, scale, n1, n2

  if n_params() lt 4 then begin
      print, 'Syntax- flux = shapelet_flux( coefs, scale, n1, n2 )'
      return, 0
  endif

  m = n_elements(coefs)

  sum = 0d
  for i = 0, m - 1 do begin
      
      if (n1[i] mod 2 eq 0) and (n2[i] mod 2 eq 0) then sum = sum + $
        2.0^((2.0 - n1[i] - n2[i]) / 2.0) * sqrt( factorial(n1[i]) * factorial(n2[i]) ) / $
        ( factorial(n1[i] / 2.0) * factorial(n2[i] / 2.0) ) * coefs[i]
      
  endfor

  flux = sqrt(!pi) * scale * sum

  return, flux
end
