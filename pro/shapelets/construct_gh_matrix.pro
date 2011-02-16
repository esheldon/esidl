;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Function to construct the Gauss-Hermite (shapelet) basis. The rows of
; this matrix contain the basis functions, e.g., GH[*,j] contains the
; jth shapelet mode sampled at the NX data points. Note that NX and X0
; are 2-element vectors.
;
; Author : Brandon C. Kelly, Steward Obs., Feb 2005
;
; Inputs :
;    NX - The dimensions of the shapelet basis array. The shapelet
;         basis images will be NX[0] x NX[1] arrays.
;    X0 - The centroid of the shapelets, a 2-element vector.
;    SCALE - The scale of the shapelets.
;    NMAX - The maximum order of the shapelets. All shapelet states
;           for which n1 + n2 =< NMAX will be generated.
;
; Output :
;    The shapelet basis, an [NX[0] * NX[1], M] array, where M is the
;    number of shapelet states. 2-d arrays of the shapelet functions
;    can be generated from the rows of GH as reform(GH[*,j], NX[0],
;    NX[1]).
;
; Optional Outputs :
;   N1, N2 - M-element array containing the order of the
;            shapelets. GH[*,j] contains the shapelet of order
;            (N1[j],N2[j]).
;
; Functions Called :
;   PHIN
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function construct_gh_matrix, nx, x0, scale, nmax, n1=n1, n2=n2

  if n_params() lt 4 then begin
      print, 'Syntax- GH = construct_gh_matrix( nx, x0, scale, nmax, n1=n1, n2=n2 )'
      return, 0
  endif

  m = long((nmax + 1) * (nmax + 2) / 2) ;number of shapelet states

  n1 = intarr(m)
  n2 = n1

  i = 0L
  for n = 0, nmax do begin      ;creat the index vectors, n1 and n2
      for n1i = 0, n do begin
          n1[i] = n1i
          n2[i] = n - n1i
          i = i + 1
      endfor
  endfor

  x1 = findgen(nx[0]) - x0[0]   ;the absissca in the horizontal direction
  x2 = findgen(nx[1]) - x0[1]   ;the absissca in the vertical direction

  GH = dblarr(nx[0]*nx[1], m)   ;the shapelet basis
;get the basis vectors and store them in GH
  for j = 0, m - 1 do GH[*,j] = phin( n1[j], n2[j], x1, x2, [scale,scale] )

  return, GH
end
