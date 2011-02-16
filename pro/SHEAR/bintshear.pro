PRO bintshear, e1, e2, uncert, $
               x, y, R, $
               rmin, rmax, binsize, $
               etan, erad, etanerr, eraderr, meanr, $
               npair=npair, $
               keepsum=keepsum, $
               etansum=etansum, eradsum=eradsum, $
               etanerrsum=etanerrsum, eraderrsum=eraderrsum, $
               wsum=wsum, npsum=npsum, rsum=rsum

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
     print,'-Syntax: bintshear, e1, e2, uncert, '
     print,'    x, y, R,'
     print,'    rmin, rmax, binsize'
     print,'    [ etan, erad, etanerr, eraderr, meanr, npair ]'
     print,''
     print,'Use doc_library,"bintshear"  for more help.'  
     return
  ENDIF 

  IF NOT keyword_set(keepsum) THEN keepsum = 0 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  IF NOT keyword_set(keepsum) THEN keepsum = 0

                                ; Default intrinsic shape noise
  vint = .32^2
                                ; Default value for arrays
  default = 1.e10

  diff = rmax - rmin
  nbin = long( diff/binsize )
  binarr, R, binsize, rmin, rmax, ind

  ;; tangential and radial in 2*theta space
  IF NOT keepsum THEN BEGIN 
      etan    = replicate(default, nbin)
      erad    = etan
      etanerr = etan
      eraderr = etan
      meanr   = etan
  ENDIF ELSE BEGIN 
      IF n_elements(npsum) EQ 0  THEN BEGIN
          etansum    = replicate(0., nbin)
          eradsum    = etansum
          etanerrsum = etansum
          eraderrsum = etansum
          wsum       = etansum
          rsum       = etansum
          npsum      = lonarr(nbin)
      ENDIF ELSE IF total(npsum) EQ 0 THEN BEGIN
          etansum    = replicate(0., nbin)
          eradsum    = etansum
          etanerrsum = etansum
          eraderrsum = etansum
          wsum       = etansum
          rsum       = etansum
          npsum      = lonarr(nbin)
      ENDIF 
  ENDELSE  
  npair = lonarr(nbin)

  FOR i=0L, nbin-1 DO BEGIN 
      IF ind(i) NE ind(i+1) THEN BEGIN 
          w = ind( ind(i):ind(i+1)-1 )
          npair[i] = n_elements(w)
          ;; Tangential and Radial in 2*theta space
          ;; Using x's and y's faster than using trig functions
          e1prime  = -( e1[w]*(x[w]^2 - y[w]^2) + e2[w]*2.*x[w]*y[w] )/R[w]^2
          e2prime  =  ( e1[w]*2.*x[w]*y[w]  -  e2[w]*(x[w]^2 - y[w]^2) )/R[w]^2

          IF NOT keepsum THEN BEGIN 
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
          ENDIF ELSE BEGIN 
              npsum[i] = npsum[i] + npair[i]
              weights = 1./( vint + uncert[w]^2 )
              
              rsum[i] = rsum[i] + total(R[w])
              etansum[i] = etansum[i] + total(e1prime*weights)
              eradsum[i] = eradsum[i] + total(e2prime*weights)
              ;; assuming <x> ~ 0.
              etanerrsum[i]  = etanerrsum[i] + total(weights^2*e1prime^2)
              eraderrsum[i]  = eraderrsum[i] + total(weights^2*e2prime^2)
              wsum[i] = wsum[i] + total(weights)
          ENDELSE 
      ENDIF 
  ENDFOR 

  return 
END 
