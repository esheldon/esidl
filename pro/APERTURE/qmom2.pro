PRO qmom2, e1, e2, xabs, yabs, cenx, ceny, Rmax, sigma, qMap, qrad, qMaperr, qraderr, nobj=nobj

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;    QMOM2
;
; PURPOSE:
;    This one's for finding <Map^2>
;
; CALL:
;
; PROCEDURE: 
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: qmom2, e1, e2, xabs, yabs, cenx, ceny, Rmax, sigma, qMap, qrad, qMaperr, qraderr, nobj=nobj'
      return
  ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up relative positons and Q
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  x = xabs - cenx
  y = yabs - ceny

  R = sqrt(x^2 + y^2)
  R2 = R^2

  wh=where( R LE Rmax, nwh)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up the weights
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  gd=where(sigma[wh] GT 0.0, nobj)
  gd = wh[gd]
  w = ( 1./sigma[gd]^2 )
  wtot = total(w)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Define Q (see Schneider 98 (97) )
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  Area = !pi*Rmax^2
  L = 1.                        ; 1. seems as good as any
  G = (1+L)*(2+L)^2/(1+2*L)/(3+2*L)

  Q = ( (1+L)*(2+L)/Area )*(R2[gd]/Rmax^2)*(1. - R2[gd]/Rmax^2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Tangential and Radial Shear
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  etan  = -( e1[gd]*(x[gd]^2 - y[gd]^2) + e2[gd]*2.*x[gd]*y[gd] )/R2[gd]
  erad  =  ( e1[gd]*2.*x[gd]*y[gd]  -  e2[gd]*(x[gd]^2 - y[gd]^2) )/R2[gd]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Weighted mean
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  qMap = 0.                     ;Note qMap is really <M^2>
  qrad = 0.

  ftan = w^2*Q*etan
  frad = w^2*Q*erad

  FOR i=0L, nobj-1 DO BEGIN 
      FOR j=0L, nobj-1 DO BEGIN
          
          qMap = qMap + (i NE j)*( ftan[i] )*( ftan[j] )
          qrad = qrad + (i NE j)*( frad[i] )*( ftan[j] )

      ENDFOR  
  ENDFOR 

  qMap = Area^2*qMap/wtot^2
  qrad = Area^2*qrad/wtot^2


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Uncertainty in the mean.  This is an estimate from Schneider 98 (97)
; converted to include estimate errors.
;
; Assumes Kurtosis is small relative to <M^2)^2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  qMaperr = sqrt(2.)*( G/wtot + .01 )
  qraderr = sqrt(2.)*( G/wtot + .01 ) ; Same order of mag as qMap
                                       ; Thats why we need lots of fields.



return
END
