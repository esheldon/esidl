PRO qmom3, e1, e2, xabs, yabs, cenx, ceny, Rmax, sigma, qMap, qrad, qMaperr, qraderr, nobj=nobj

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

  wh=where( R LE Rmax, nobj)
  IF nobj EQ 0 THEN return

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up the weights
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  gd=where(sigma[wh] GT 0.0, nobj)
  nobj=long(nobj)
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

  ftan = w^2*Q*etan*sqrt(2.*Area^2/wtot^2)
  frad = w^2*Q*erad*sqrt(2.*Area^2/wtot^2)


  IF (nobj EQ 1) THEN BEGIN

      qMap = 0.
      qrad = 0.

  ENDIF ELSE BEGIN 

      nind = nobj*(nobj-1L)/2L

      aa = indgen(nobj)
      ii=intarr(nind)
      jj=ii

      s1=-1L
      s2=-1L
      FOR i=0L, nobj-2 DO BEGIN
              
          s1 = s2 + 1L
          s2 = s1 + nobj - (i+2)
          ii[ s1:s2 ] = i
          jj[ s1:s2 ] = aa[ i+1:nobj-1 ]
      ENDFOR 

      step = 100000L
      IF nind GT step THEN BEGIN 

          ;; If too big, Break into chunks
          nstep = long(nind/step)
          left = nind MOD step

          aa = lindgen(step)
          FOR istep = 0, nstep-1 DO BEGIN 
          
              aa = aa+istep*step
              useii = ii[aa]
              usejj = jj[aa]
          
              fsend = ftan[useii]*ftan[usejj]
              qMap = qMap + total( temporary(fsend) )
              fsend=0
              fsend = frad[useii]*frad[usejj]
              qrad = qrad + total( temporary(fsend) )

          ENDFOR 
          IF left NE 0 THEN BEGIN 
              aa = aa + istep*step ;Using exit value of istep
              useii = ii[aa[0:left-1]]
              usejj = jj[aa[0:left-1]]
              fsend = ftan[useii]*ftan[usejj]
              qMap = qMap + total( temporary(fsend) )
              fsend=0
              fsend = frad[useii]*frad[usejj]
              qrad = qrad + total( temporary(fsend) )
          ENDIF 
      ENDIF ELSE BEGIN 
          fsend = ftan[ii]*ftan[jj]
          qMap = qMap + total( temporary(fsend) )
          fsend=0
          fsend = frad[ii]*frad[jj]
          qrad = qrad + total( temporary(fsend) )
      ENDELSE 
  ENDELSE 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Uncertainty in the mean.  This is an estimate from Schneider 98 (97)
; converted to include estimate errors.
;
; Assumes Kurtosis is small relative to <M^2)^2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  qMaperr = sqrt(2.)*( G/wtot + .01 )
  qraderr = sqrt(2.)*( G/wtot + .01 ) ; Same order of mag as qMap
                                       ; Thats why we need lots of fields.

  ftan=0
  frad=0
  fsend=0
  wh=0
  g=0
  ii=0
  jj=0                          ;Free the memory

return
END
