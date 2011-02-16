PRO qmom, e1, e2, xabs, yabs, cenx, ceny, Rmax, sigma, qtan, qrad, qtanerr, qraderr, qtansig, qradsig, nobj=nobj

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:
;    QMOM
;
; PURPOSE:
;    
;
; CALL:
;
; PROCEDURE: 
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: qmom, e1, e2, xabs, yabs, cenx, ceny, Rmax, sigma, qtan, qrad, qtanerr, qraderr, qtansig, qradsig, nobj=nobj'
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

  g=where(sigma[wh] GT 0.0, nobj)
  g = wh[g]
  w = ( 1./sigma[g]^2 )
  wtot = total(w)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Define Q
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  Q = (6./!pi)*R2[g]*(Rmax^2 - R2[g])/Rmax^6

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Tangential and Radial Shear
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  etan  = -( e1[g]*(x[g]^2 - y[g]^2) + e2[g]*2.*x[g]*y[g] )/R2[g]
  erad  =  ( e1[g]*2.*x[g]*y[g]  -  e2[g]*(x[g]^2 - y[g]^2) )/R2[g]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Weighted mean
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  qtan = total( w*Q*etan )/wtot
  qrad = total( w*Q*erad )/wtot

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Uncertainty in the mean.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  qtanerr = sqrt( total( w^2*(Q*etan - qtan)^2 )/wtot^2 )
  qraderr = sqrt( total( w^2*(Q*erad - qrad)^2 )/wtot^2 )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Variance about the mean.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  qtansig = sqrt( total( w*(Q*etan - qtan)^2)/wtot )
  qradsig = sqrt( total( w*(Q*erad - qrad)^2)/wtot )



return
END
