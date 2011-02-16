PRO mygcirc_survey,lam1,eta1,lam2,eta2,dis,theta,radians_out=radians_out
;+
; NAME:
;     MYGCIRC_SURVEY
; PURPOSE:
;     Computes rigorous great circle arc distances.  
;     Adapted from gcirc in idlastron library to work with SDSS survey
;     coordinates.
;
; PROCEDURE CALLS:
;      ISARRAY()
;
; NOTES:
;       (1) If LAM1,ETA1 are scalars, and RA2,ETA2 are vectors, then DIS is a
;       vector giving the distance of each element of LAM2,ETA2 to LAM1,ETA1.
;       Similarly, if LAM1,ETA1 are vectors, and LAM2, ETA2 are scalars, then DIS
;       is a vector giving the distance of each element of LAM1, ETA1 to 
;       LAM2, ETA2.    If both LAM1,ETA1 and LAM2,ETA2 are vectors then DIS is a
;       vector giving the distance of each element of LAM1,ETA1 to the 
;       corresponding element of LAM2,ETA2.    If the input vectors are not the 
;       same length, then excess elements of the longer ones will be ignored.
;
;       (2) Coordinates closer together than a few milliarcsec cannot
;       be distinguished.  If you are in this realm, you should be
;       using special-purpose algorithms.
;
;
;   MODIFICATION HISTORY:
;      Written in Fortran by R. Hill -- SASC Technologies -- January 3, 1986
;      Translated from FORTRAN to IDL, RSH, STX, 2/6/87
;      Vector arguments allowed    W. Landsman    April 1989
;      Prints result if last argument not given.  RSH, RSTX, 3 Apr. 1998
;      Converted to IDL V5.0                      April 1998
;      Remove ISARRAY(), V5.3 version        W. Landsman   August 2000
;      Stripped down, only degrees allowed. E.S.S.
;      Added theta calculation. E.S.S. UofMich 2002
;-
  On_error,2                    ;Return to caller

  npar = N_params()
  IF npar LT 5 THEN BEGIN
      print,'-Syntax: mygcirc_survey,lam1,eta1,lam2,eta2,dis,theta,radians_out=radians_out'
      RETURN
  ENDIF
  
  scalar = (size(lam1,/N_Dimen) EQ 0) and (size(lam2,/N_dimen) EQ 0)
  IF scalar THEN BEGIN
      IF (lam1 eq lam2) and (eta1 eq eta2) THEN BEGIN
          dis = 0.0d0
          theta = 0.0d0
          print,'Positions are equal:  ', lam1, eta1
          return
      ENDIF
  ENDIF
  
  d2r    = !DPI/180.0d0
  r2d    = 180d/!DPI
  

;  sinlam1 = sin(lam1*d2r)
;  coslam1 = cos(lam1*d2r)

;  sinlam2 = sin(lam2*d2r)
;  coslam2 = cos(lam2*d2r)

  lam1rad = lam1*d2r
  sinlam1 = sin(lam1rad)
  coslam1 = cos(lam1rad)
  lam1rad = 0

  lam2rad = lam2*d2r
  sinlam2 = sin(lam2rad)
  coslam2 = cos(lam2rad)
  lam2rad = 0


  etadiff = (eta2-eta1)*d2r
; etadiff = (eta1-eta2)*d2r
  cosetadiff = cos(etadiff)

  cosdis = sinlam1*sinlam2 + coslam1*coslam2*cosetadiff
  ww=where(cosdis LT -1., nww)
  IF nww NE 0 THEN cosdis[ww] = -1.
  ww=where(cosdis GT 1., nww)
  IF nww NE 0 THEN cosdis[ww] = 1.
 
  dis = acos(cosdis)
  IF NOT keyword_set(radians_out) THEN dis = dis*r2d

  IF arg_present(theta) THEN BEGIN ;; calculate theta
      IF n_elements(lam1) NE 1 THEN BEGIN 
          message,'Can only find theta for # in lam1,eta1 array = 1'
      ENDIF 
      theta = atan( sin(etadiff), $
                    (sinlam1*cosetadiff - coslam1*sinlam2/coslam2) ) - !dpi/2.
  ENDIF 
  
  
  RETURN
 END

