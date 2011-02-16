PRO sum_shear_fix2, e1, e2, uncert, $
                   R, x, y, $
                   rmin, rmax, binsize, $
                   etansum, eradsum, $
                   etanerrsum, eraderrsum, $
                   wsum, npsum, rsum, Sshsum, npair

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
     print,'-Syntax: sum_shear, e1, e2, uncert, '
     print,'    xabs, yabs, cenx, ceny,'
     print,'    rmin, rmax, binsize'
     print,'    [ etansum, eradsum, etanerrsum, eraderrsum, '
     print,'      wsum, npsum, rsum, npair ]'
     print,''
     print,'Use doc_library,"bintshear"  for more help.'  
     return
  ENDIF 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
                                
  vint = .32^2                  ; Default intrinsic shape noise

  IF n_elements(R) EQ 1 THEN R=[R]

  diff = rmax - rmin
  nbin = long( diff/binsize )
  binarr, R, binsize, rmin, rmax, ind

  ;; tangential and radial in 2*theta space

  npair = lonarr(nbin)
  IF n_elements(npsum) EQ 0  THEN BEGIN
      etansum    = fltarr(nbin)
      eradsum    = etansum
      etanerrsum = etansum
      eraderrsum = etansum
      wsum       = etansum
      rsum       = etansum
      npsum      = etansum
      Sshsum     = 0.           ;Scalar
  ENDIF 

  FOR i=0L, nbin-1 DO BEGIN 
      IF ind(i) NE ind(i+1) THEN BEGIN ; Make sure bin is not empty

          w = ind( ind(i):ind(i+1)-1 )
          npair[i] = n_elements(w)

          ;; Tangential and Radial in 2*theta space
          ;; Using x's and y's faster than using trig functions
          e1prime  = -( e1[w]*(x[w]^2 - y[w]^2) + e2[w]*2.*x[w]*y[w] )/R[w]^2
          e2prime  =  ( e1[w]*2.*x[w]*y[w]  -  e2[w]*(x[w]^2 - y[w]^2) )/R[w]^2

          npsum[i] = npsum[i] + npair[i]
          weights = 1./( vint + uncert[w]^2 )
              
          rsum[i] = rsum[i] + total(R[w])
          etansum[i] = etansum[i] + total(e1prime*weights)
          eradsum[i] = eradsum[i] + total(e2prime*weights)
          ;; assuming <x> ~ 0.
          etanerrsum[i]  = etanerrsum[i] + total(weights^2*e1prime^2)
          eraderrsum[i]  = eraderrsum[i] + total(weights^2*e2prime^2)

          Sshsum = Sshsum + total( weights*(1.-vint*weights*e1prime^2) )

          wsum[i] = wsum[i] + total(weights)

      ENDIF 
  ENDFOR 

  return 
END 
