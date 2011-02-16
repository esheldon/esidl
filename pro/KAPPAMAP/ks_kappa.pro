PRO ks_kappa, e1, e2, uncert, $
              xabs, yabs, cenx, ceny, $
              rmax, slength, $
              kappa, radial, wsum, npair, kerr1, kerr2, rerr1, rerr2, err3, $
              fgal=fgal, dens=dens

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
     print,'-Syntax: ks_kappa, e1, e2, uncert, xabs, yabs, cenx, ceny, rmax, slength, kappa, radial, wsum, npair'
     print,''
     print,'Use doc_library,"ks_kappa"  for more help.'  
     return
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  x = xabs - cenx            ;Relative positions
  y = yabs - ceny

  R=sqrt( x^2 + y^2 )
  w = where( R LE rmax, npair)  ; Find stuff within rmax
;  w = where( R LE rmax AND R GT 10./3600., npair)

  IF npair EQ 0 THEN return

  R2 = R[w]^2
  ratio = -R2/2./slength^2
  xy = x[w]*y[w]
  xminy = x[w]^2 - y[w]^2

  Q = ( 1. - (1. - ratio )*exp(ratio) )/R2/!pi
;  Q = (6./!pi)*R2*(rmax^2 - R2)/rmax^6
  Q2 = Q^2

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Tangential and Radial in 2*theta space
  ;; Using x's and y's faster than using trig
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  etan  = -( e1[w]*xminy + e2[w]*2.*xy )/R2
  erad  =  ( e1[w]*2.*xy - e2[w]*xminy )/R2
  
  vint = .32^2
  weight = 1./(vint + uncert[w]^2)
  weight2 = weight^2
  
  kappa = total(Q*etan*weight)
  radial= total(Q*erad*weight)

  ;; For calculating the error bars
  wsum  = total(weight)

  kerr1 = total(Q2*etan^2*weight2)
  kerr2 = total(Q*etan*weight2)

  rerr1 = total(Q2*erad^2*weight2)
  rerr2 = total(Q*erad*weight2)

  err3  = total(weight2)

  IF n_elements(fgal) NE 0 THEN BEGIN 
      
      dlength = slength
      xf = fgal.dec - cenx
      yf = fgal.ra  - ceny

      RL = sqrt(xf^2 + yf^2)
      
      wL = where( RL LE rmax, nwL)
      IF nwL NE 0 THEN BEGIN 
          dens = total( exp(-RL^2/2./dlength^2) )
      ENDIF ELSE dens = 0.
  ENDIF 

  return 
END 
