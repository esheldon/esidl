FUNCTION sdss_coords::init
  return 1
END 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
+
;
; NAME:
;    SDSS_COORDS::EQ2CSURVEY
;       
; PURPOSE:
;    Convert from ra, dec to the corrected clambda, ceta 
;    SDSS survey coordinate system.  It is corrected so that the
;    longitude eta ranges from [-180.0, 180.0] and the latitude
;    lambda ranges from [-90.0,90.0].  The standard lambda/eta 
;    both range from [-180.0,180.0] which doesn't make sense.
;    NOTE: lambda is often referred to as longitude but this
;    is incorrect since it has poles at [-90,90]
;
; CALLING SEQUENCE:
;    sc->eq2csurvey, ra, dec, clambda, ceta
;
; INPUTS: 
;    ra: Equatorial latitude in degrees 
;    dec: Equatorial longitude in degrees
;
; OPTIONAL INPUTS:
;    None
;
; KEYWORD PARAMETERS:
;    None
;       
; OUTPUTS: 
;    clambda: Corrected Survey longitude (actually lattitude) in degrees
;    ceta: Corrected Survey latitude (actually logitude) in degrees
;
; OPTIONAL OUTPUTS:
;    None
;
; CALLED ROUTINES:
;    ATBOUND
;    ATBOUND2
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    Written: 26-Sep-2002  Erin Scott Sheldon
;       
;                                      
;docend::sdss_coords::eq2csurvey
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO eq2csurvey, ra_in, dec_in, clambda, ceta

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: eq2csurvey, ra, dec, clambda, ceta'
      print,' ra, dec in degrees'
      return
  ENDIF 
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Some parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  surveyCenterRa  =  185.0
  surveyCenterDec =   32.5

  deg2Rad = !dpi/180.
  rad2Deg = 180./!dpi
 
  node = surveyCenterRa - 90.
  node = node * deg2Rad

  etaPole = surveyCenterDec * deg2Rad

  nra = n_elements(ra_in)
  ndec = n_elements(dec_in)
  IF nra NE ndec THEN BEGIN 
      print,'ra and dec must be same size'
      return
  ENDIF 

  w=where((ra_in GT 360d) OR (ra_in LT 0d), nw)
  IF nw NE 0 THEN message,'RA must be within [0,360]'
  w=where((dec_in GT 90d) OR (dec_in LT -90d),nw)
  IF nw NE 0 THEN message,'DEC must be within [-90,90]'

  ;; Convert to radians
  ra = double(ra_in*deg2Rad)
  dec = double(dec_in*deg2Rad)

  x1 = cos(ra-node)*cos(dec)
  y1 = sin(ra-node)*cos(dec)
  z1 = sin(dec)

  ;; free memory
  ra=0
  dec=0

  ;; calculate corrected survey coordinates

  clambda = -asin( temporary(x1) )
  ceta = atan( temporary(z1), temporary(y1) ) - etaPole

  ;; convert to degrees
  clambda = clambda * rad2Deg
  ceta = ceta * rad2Deg

  ;; make sure ceta is between -180.0 and 180.0
  self->longbound, ceta, -180.0, 180.0

  return

END



PRO sdss_coords::longbound, longitude, min, max

  w=where(longitude LT min,nw)
  WHILE nw NE 0 DO BEGIN 
      longitude[w] = longitude[w] + 360.0
      w=where(longitude LT min,nw)
  ENDWHILE 

  w=where(longitude GE max,nw)
  WHILE nw NE 0 DO BEGIN 
      longitude[w] = longitude[w] - 360.0
      w=where(longitude GT max,nw)
  ENDWHILE 

  return
END 


PRO sdss_coords__define

  struct = {$
             sdss_coords, $
             sdss_util_dummy: 0 $
           }

END 
