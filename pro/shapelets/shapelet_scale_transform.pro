;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Function to compute the shapelet coefficients if some new scale were
; to be used, from the current coefficients and shapelet scale. See
; the appendix of Refregier (MNRAS, 2003) for more detail.
;
; Author : Brandon C. Kelly, Steward Obs., Feb 2005
;
; Inputs :
;   COEFS - The shapelet coefficients for shapelet of scale SCALE1.
;   SCALE1 - The scale corresponding to COEFS
;   SCALE2 - The new scale that one wants the shapelet coefficients
;            for.
;   N10, N20 - The N1 and N2 arrays that contain the order of the
;              shapelet coefficients. More specifically,
;              (N10[j],N20[j]) is the order of the jth shapelet. These
;              vectors should be the same size as COEFS.
; Output :
;   The shapelet coefficients corresponding to shapelets of scale
;   SCALE2.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function shapelet_scale_transform, coefs, scale1, scale2, n10, n20

  if n_params() lt 5 then begin
      print, 'Syntax- new_coefs = shapelet_scale_transform( coefs, scale1, '
      print, '                                              scale2, n10, n20 )'
      return, 0
  endif

  nmax = max(n10 + n20)         ;get max shapelet order

  b1 = (scale1^2 - scale2^2) / (scale1^2 + scale2^2)
  b2 = 2.0 * scale1 * scale2 / (scale1^2 + scale2^2)

  T = dblarr(nmax+1, nmax+1)    ;the 1-d transformation matrix

;create the transformation matrices
  for n1 = 0, nmax do begin
      
      for n2 = 0, nmax do begin
                                ;n1 and n2 must be both odd or even,
          IF ((n1 MOD 2) EQ (n2 MOD 2)) THEN BEGIN
              
              case 1 of         ;set up l vector
                  min([n1, n2]) eq 0 : l = 0
                  min([n1, n2]) eq 1 : l = 1
                  (n1 mod 2) eq 0 : l = 2 * indgen( min([n1, n2]) / 2 + 1 ) ;l must be even
                  (n1 mod 2) eq 1 : l = 2 * indgen( ceil( min([n1, n2]) / 2.0 ) ) + 1 ;l must be odd
              endcase
              
              entry = (-1.0)^( (n2 - l) / 2.0 ) * sqrt( factorial(n1) * factorial(n2) ) / $
                ( factorial( (n1 - l) / 2.0 ) * factorial( (n2 - l) / 2.0 ) * factorial(l) ) * $
                (b1 / 2.0)^((n1 + n2) / 2.0 - l) * b2^(l + 0.5)
              
              T[n2,n1] = total(entry, /nan) ;n2 indexes the column, n1 the row

          ENDIF
          
      endfor
      
  endfor

;get transformation matrix for 2-d from the 1-d matrices
  T1 = T[*,n10]
  T1 = T1[n10,*]
  T2 = T[*,n20]
  T2 = T2[n20,*]

  T = T1 * T2

;now transform the coefficients to new scale
  new_coefs = T ## transpose(coefs)

  return, reform(new_coefs)
end
