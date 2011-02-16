PRO mygcirc,ra1,dc1,ra2,dc2,dis,theta,radians_out=radians_out
;+
; NAME:
;     MYGCIRC
; PURPOSE:
;     Computes great circle distances and angles between sets of points.
; INPUTS:
;   ra1,dec1: First list of positions in degrees.
;   ra2, dec2: Second list of positions in degrees.
;       The lengths of these sets of coordss must make sense for calculations 
;       such as as sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos(ra1-ra2)
;       A typical use is a single ra1,dec1 and a longer list of ra2,dec2.
; KEYWORDS:
;   /radians_out: Return dis in radians instead of degrees.
; OUTPUTS:
;   distance: The distance in degrees unless /radians_out is set
;   theta: The angle between the second point and the (XXX fill this in).
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
 On_error,2                            ;Return to caller

 npar = N_params()
 IF npar LT 5 THEN BEGIN
   print,'-Syntax: mygcirc,ra1,dc1,ra2,dc2,dis,theta,radians_out=radians_out'
   RETURN
 ENDIF

 scalar = (size(ra1,/N_Dimen) EQ 0) and (size(ra2,/N_dimen) EQ 0)
 IF scalar THEN BEGIN
    IF (ra1 eq ra2) and (dc1 eq dc2) THEN BEGIN
       dis = 0.0d0
       theta = 0.0d0
       print,'Positions are equal:  ', ra1, dc1
       return
    ENDIF
 ENDIF

 d2r    = !DPI/180.0d0
 r2d    = 180.0d0/!DPI

 sindc1 = sin(dc1*d2r)
 cosdc1 = cos(dc1*d2r)

 sindc2 = sin(dc2*d2r)
 cosdc2 = cos(dc2*d2r)

; radiff = (ra2-ra1)*d2r
 radiff = (ra1-ra2)*d2r
 cosradiff = cos(radiff)

 cosdis = sindc1*sindc2 + cosdc1*cosdc2*cosradiff
 ww=where(cosdis LT -1., nww)
 IF nww NE 0 THEN cosdis[ww] = -1.
 ww=where(cosdis GT 1., nww)
 IF nww NE 0 THEN cosdis[ww] = 1.

 dis    = acos(cosdis)
 IF NOT keyword_set(radians_out) THEN dis = dis*r2d

 IF arg_present(theta) THEN BEGIN ;; calculate theta
     IF n_elements(ra1) NE 1 THEN BEGIN 
         message,'Can only find theta for # in ra1,dec1 array = 1'
         return
     ENDIF 
     theta = atan( sin(radiff), $
                   (sindc1*cosradiff - cosdc1*sindc2/cosdc2) ) - !dpi
 ENDIF 

 RETURN
 END

