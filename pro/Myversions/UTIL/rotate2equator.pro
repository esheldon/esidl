
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    ROTATE2EQUATOR
;       
; PURPOSE:
;    Rotate coord. system so that the input SDSS great circle inclination
;    is centered on the equator, where the coordinate system is rectangular 
;    within the width of a stripe (2.5 degrees)
;
; CALLING SEQUENCE:
;    rotate2equator, ra, dec, inc, ra_prime, dec_prime, plot=plot
;
; INPUTS: 
;    ra,dec: ra and dec vectors. Must be same size.
;    inc: the inclination angle between the SDSS great circle and
;         the equator (after a rotation of 5 degrees in ra)
;
; OPTIONAL INPUTS:
;    None
;
; KEYWORD PARAMETERS:
;    /plot: make plot of before and after ra,dec vectors
;    /noresetdomain: allow the new ra to go -180,180 instead
;            of resetting it to [0,360]
;       
; OUTPUTS: 
;    ra_prime,dec_prime: vectors in new coordinate system.
;
; OPTIONAL OUTPUTS:
;    None
;
; CALLED ROUTINES:
;    
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    10-Mar-2001 Erin Scott Sheldon UofMich
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO rotate2equator, ra, dec, inc, ra_prime, dec_prime, plot=plot, noresetdomain=noresetdomain, rotsouth=rotsouth

  IF N_params() LT 3 THEN BEGIN 
     print,'-Syntax: rotate2equator, ra, dec, inc, ra_prime, dec_prime, plot=plot, noresetdomain=noresetdomain'
     print,''
     print,'Use doc_library,"rotate2equator"  for more help.'  
     return
  ENDIF 

  nra =  n_elements(ra)
  ndec = n_elements(dec)

  IF nra NE ndec THEN message,'ra and dec must be same size'

  r2d = 180d/!dpi
  d2r = !dpi/180d

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Survey coordinates are centered on 185.0, so rotate around 
  ;; z-axis by 5.0, then we can rotate around the y-axis by
  ;; the inclination angle, which is defined relative to 
  ;; survey coordinates.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  deltaz = 5.0*d2r              
  deltay = inc*d2r               ;Rotate around y axis by the inclination

  R = 1d
  theta = (90d - dec)*d2r
  phi = ra*d2r

  rho = R*sin(theta)
  
  x = rho*cos(phi)
  y = rho*sin(phi)
  z = R*cos(theta)

  xprime=dblarr(nra)
  yprime=xprime
  zprime=yprime

  ;; y rotation matrix

         ;col1      col2    col3
  yrot = $
    [ [ cos(deltay), 0d, -sin(deltay)], $ ;row1
      [ 0d,          1d,           0d], $ ;row2
      [ sin(deltay), 0d,  cos(deltay)]  ] ;row3


  ;; z rotation matrix

         ;col1      col2    col3
  zrot = $
    [ [  cos(deltaz), sin(deltaz), 0d], $ ;row1
      [ -sin(deltaz), cos(deltaz), 0d], $ ;row2
      [      0d,          0d,      1d]  ] ;row3

  ;; z rotation to recenter southern stripes
  ;; as if they were northern, so don't cross 0/360 mark

  deltaz2 = 180.*d2r
  zrot2 = $
    [ [  cos(deltaz2), sin(deltaz2), 0d], $ ;row1
      [ -sin(deltaz2), cos(deltaz2), 0d], $ ;row2
      [      0d,          0d,        1d]  ] ;row3

  ;; rotate
  
  IF keyword_set(rotsouth) THEN BEGIN 
      FOR i=0L, nra-1 DO BEGIN 
          
          vector = [x[i], y[i], z[i]]
          vector_prime = reform( zrot2##(yrot##(zrot##transpose(vector))) )
          
          xprime[i] = vector_prime[0]
          yprime[i] = vector_prime[1]
          zprime[i] = vector_prime[2]
          
      ENDFOR 
  ENDIF ELSE BEGIN 
      FOR i=0L, nra-1 DO BEGIN 
          
          vector = [x[i], y[i], z[i]]
          vector_prime = reform( yrot##(zrot##transpose(vector)) )
          
          xprime[i] = vector_prime[0]
          yprime[i] = vector_prime[1]
          zprime[i] = vector_prime[2]
          
      ENDFOR 
  ENDELSE 

  ;; now get angles from x,y,z
  phi_prime = atan( yprime, xprime )
  theta_prime = acos( zprime/R )

  dec_prime = 90d - theta_prime*r2d
  ra_prime = phi_prime*r2d

  ;; put back in the 5 degrees by which we rotated
  ra_prime = ra_prime + 5d

  ;; Make sure ra in [0,360] unless requested otherwise
  IF NOT keyword_set(noresetdomain) THEN BEGIN 
      w=where(ra_prime LT 0.0,nw)
      IF nw NE 0 THEN ra_prime[w] = ra_prime[w]+360.
      w=where(ra_prime GT 360.,nw)
      IF nw NE 0 THEN ra_prime[w] = ra_prime[w]-360.
  ENDIF 

  IF keyword_set(plot) THEN BEGIN 
      !p.multi=[0,0,2]
      
      plot,ra,dec,psym=3,xtitle='ra',ytitle='dec'
      plot, ra_prime, dec_prime, psym=3,xtitle="ra'",ytitle="dec'"

      !p.multi=0
  ENDIF 

  IF nra EQ 1 THEN BEGIN 
      ra_prime=ra_prime[0]
      dec_prime=dec_prime[0]
  ENDIF 

END 
