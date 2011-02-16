PRO tan_shear, e1, e2, uncert, $
               xabs, yabs, cenx, ceny, $
               rmin, rmax, binsize, $
               etan, erad, etanerr, eraderr, meanr, npair

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;       
; PURPOSE:
;	
;
; CALLING SEQUENCE:
;      
;                 
;
; INPUTS: 
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
;       
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; CALLED ROUTINES:
; 
; PROCEDURE: 
;	
;	
;
; REVISION HISTORY:
;	
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() lt 9 THEN BEGIN 
     print,'-Syntax: tan_shear, e1, e2, uncert, '
     print,'    xabs, yabs,'
     print,'    rmin, rmax, binsize'
     print,'    [ etan, erad, etanerr, eraderr, meanr, npair ]'
     print,''
     print,'Use doc_library,"bintshear"  for more help.'  
     return
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                ; Default intrinsic shape noise
  vint = .32^2
                                ; Default value for arrays
  default = 1.e10

  x = xabs - cenx
  y = yabs - ceny

  R=sqrt( x^2 + y^2 ) > double(1.e-6)

  diff = rmax - rmin
  nbin = long( diff/binsize )
  binarr, R, binsize, rmin, rmax, ind

  ;; tangential and radial in 2*theta space
  etan    = replicate(default, nbin)
  erad    = etan
  etanerr = etan
  eraderr = etan
  meanr   = etan

  npair = lonarr(nbin)

  FOR i=0L, nbin-1 DO BEGIN 
      IF ind(i) NE ind(i+1) THEN BEGIN 

          w = ind( ind(i):ind(i+1)-1 )
          npair[i] = n_elements(w)

          ;; Tangential and Radial in 2*theta space
          ;; Using x's and y's faster than using trig functions
          e1prime  = -( e1[w]*(x[w]^2 - y[w]^2) + e2[w]*2.*x[w]*y[w] )/R[w]^2
          e2prime  =  ( e1[w]*2.*x[w]*y[w]  -  e2[w]*(x[w]^2 - y[w]^2) )/R[w]^2

          IF npair[i] EQ 1 THEN BEGIN 
              meanr[i] = R[w]
              etan[i] = e1prime
              etanerr[i] = uncert[w]
              erad[i] = e2prime
              eraderr[i] = uncert[w]
          ENDIF ELSE BEGIN 
              IF npair[i] LT 5 THEN meanr[i] = mean( R[w] ) $
              ELSE meanr[i] = median( R[w] )

              err = sqrt( vint + uncert[w]^2 )
              wmom, e1prime, err, wmean, wsig, werr
              etan[i] = wmean
              etanerr[i] = werr

              err = sqrt( vint + uncert[w]^2 )
              wmom, e2prime, err, wmean, wsig, werr
              erad[i]    = wmean
              eraderr[i] = werr
          ENDELSE 

      ENDIF 
  ENDFOR 

  return 
END 
