FUNCTION  coord_map, x, y, $
                xprime,yprime,  $
                a=a, b=b, $
                flipx=flipx, flipy=flipy

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: coord_map
;       
; PURPOSE: find the coordinate map between primed and unprimed coordinates.
;	
;
; CALLING SEQUENCE:  coord_map, x, y, 
;                               xprime, yprime,
;                               a=a, b=b,
;                               flipx=flipx, flipy=flipy 
;
; INPUTS: x, y  A list of x and y positons in unprimed frame.  MUST BE
;               AT LEAST 3 POSITIONS IN THE LIST.  The more the merrier 
;               to help break any degeneracy.
;         xprime, yprime. a list of x,y positions in primed frame.  Must
;               correspond exactly with those from unprimed frame.
;         NOTE:  Third (or more) sets are to find stretch factors and 
;           also to break any degeneracy of the x or y positions.
;           !!!! HAVE NOT IMPLEMENTED THIS YET
;
; INPUT KEYWORD PARAMETERS:
;         a: can input the x-direction stretch factor.
;         b: same for y
;         flipx: if the primed x-axis is mirror flipped from unprimed
;         flipy: same for y 
;
;       
; OUTPUTS: map:  the mapping in a structure.  See below. 
;
; OPTIONAL OUTPUTS:
;
; CALLED ROUTINES:
; 
; PROCEDURE: 
;	Assume that the map takes the form of a translation plus a rotation
; and a stretching (or shrinking).
;
;  [ a*xp ]     [ costheta   sintheta ] [ x ]     [ c ]
;  [      ]  =  [                     ] [   ]  +  [   ]
;  [ b*yp ]     [ -sintheta  costheta ] [ y ]     [ d ]
;	
; Where xp is short for xprime
;
; - a and be are the stretching factors (the same for square pixels)
; NOTE:  a and b will be assumed to be positive numbers unless flipx or 
;      flipy are set.  This implies a mirror flip of an axis.
; - theta is the rotation angle
; - c and d are the traslations
;
; a^2, b^2, c, d, costheta and sintheta are given by the following formulae:
;
; let x12 = x1 - x2.    r12^2 = x12^2  +  y12^2 then:
;
;
;             [ r12^2   yp12^2 ]                     [ xp12^2   r12^2 ]
;         Det [                ]                 Det [                ]
;             [ r13^2   yp13^2 ]                     [ xp13^2   r13^2 ]
; a^2 =   -------------------------     b^2 =    -------------------------
;                   Det[A]                                 Det[A]
;
; 
;                              [ xp12^2   yp12^2 ]
;                          A = [                 ]
;                              [ xp13^2   yp13^2 ]
;
;
; 
;                 [ a*xp12   y12 ]                   [ x12   a*xp12 ]
;             Det [              ]               Det [              ]
;                 [ b*yp12  -x12 ]                   [ y12   b*yp12 ]
; costheta =  ----------------------- sintheta = -----------------------
;                      Det[B]                             Det[B]
;
;
;                                [ x12    y12 ]
;                            B = [            ]
;                                [ y12   -x12 ]
;
; Then multiply a or b by a minus if flipx or flipy are set.  Finally:
;
; c = a*xp1 - ( costheta*x1 + sintheta*y1 )
;
; d = b*yp1 - ( -sintheta*x1 + costheta*y1 )
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Applying the map:
;
; Use the fact that in IDL array1*array2 is not the dot product (## is)
; This actually multiplies the numbers element by element. Thus, I have
; defined arrays in the output structure map:
; 
; map.stretch = [a,b]  
; map.trans = [c,d] 
;           [ costheta  sintheta ]
; map.rot = [                    ]
;           [ -sineheta costheta ]
;
; Then, if you want to transform from unprimed to primed coords:
;
; IDL> pos = [x, y]
; IDL> result = ( map.rot##pos + map.trans )/map.stretch
;
; To transform from primed to unprimed coords:
;
; IDL> posprime = [xp, yp]
; IDL> result = invert(map.rot)##( map.stretch*posprime - map.trans)
;
;
; REVISION HISTORY: Author:  Erin Scott Sheldon  U. of Michigan 3/25/99
;	
;
; THINGS TO DO:  Implement correction of degeneracy effects.       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Define the map structure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  map = {$
    stretch: fltarr(2), $
    rot: fltarr(2,2), $
    trans: fltarr(2) }

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check for the parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if N_params() lt 4 then begin
     print,'-Syntax: coord_map, x, y, xprime, yprime, a=a, b=b, flipx=flipx, flipy=flipy'

     print,''
     print,'Use doc_library,""  for more help.'  
     return,map
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check some keywords and set some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; How small certain quantities can be.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  tol = 1.e-6

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; do we calculate the stretches?
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if ( not keyword_set(a) ) then doa = 1 
  if ( not keyword_set(b) ) then dob = 1 
  
  ;;;;; check how many entries there are 

  nx = n_elements(x)
  ny = n_elements(y)
  npx = n_elements(xprime)
  npy = n_elements(yprime)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; quit if arrays are of different size.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if (nx ne ny   or   nx ne npx   or   nx ne npy) then begin
     print,'Position arrays must be of same length.'
     return,map
  endif

  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; define some quantities
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  x1 = x[0]
  y1 = y[0]
  x2 = x[1]
  y2 = y[1]
  xp1 = xprime[0]
  yp1 = yprime[0]
  xp2 = xprime[1]
  yp2 = yprime[1]  

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; check if we are to do either stretch factor and if we have
  ; enough elements to do so
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  pos=2
  if ( nx ge 3 and (doa or dob) ) then begin
     x3 = x[pos]
     y3 = y[pos]
     xp3 = xprime[pos]
     yp3 = yprime[pos]
     pos=pos+1
  endif else if ( nx lt 3 and (doa or dob) ) then begin
       print,'Cannot find stretch factors without 3 points.'
       return,map
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; difference quantities.  Check for degeneracies on positions.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  noquit=1
  WHILE (noquit) do BEGIN
     x12 = x1 - x2
     y12 = y1 - y2
     xp12 = xp1 - xp2
     yp12 = yp1 - yp2
     if ( abs(x12) lt tol or abs(y12) lt tol or abs(xp12) lt tol $
     or abs(yp12) lt tol ) then begin

        ;;; Make sure we have enough points to break degeneracy ;;;
        if ( pos gt nx-1 ) then begin
           print,'Degeneracy in positions.  Need more points.'
           return,map
        endif 
        x2 = x[pos]
        y2 = y[pos]
        xp2 = xprime[pos]
        yp2 = yprime[pos]
        pos = pos + 1
     endif else noquit = 0
  ENDWHILE

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Define the coefficient matrix used to find costheta and sintheta
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  coeffang = fltarr(2,2)
  coeffang[0,0] = x12
  coeffang[0,1] = y12
  coeffang[1,0] = y12
  coeffang[1,1] = -x12
  detang = determ( coeffang )

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Need this matrix to find either a or b stretch factors 
  ; Check for degeneracies
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if ( doa or dob ) then begin
     coeffab = fltarr(2,2)
     noquit = 1

     WHILE (noquit) do BEGIN
     
       x13 = x1 - x3
       y13 = y1 - y3
       xp13 = xp1 - xp3
       yp13 = yp1 - yp3
       r12 = x12^2 + y12^2
       r13 = x13^2 + y13^2

       coeffab[0,0] = xp12^2
       coeffab[0,1] = xp13^2
       coeffab[1,0] = yp12^2
       coeffab[1,1] = yp13^2
       detab = coeffab[0,0]*coeffab[1,1] - coeffab[0,1]*coeffab[1,0]

       if (abs(x13) lt tol or abs(y13) lt tol or abs(xp13) lt tol $
       or abs(yp13) lt tol or abs(detab) lt tol) then begin
           if (pos gt nx-1) then begin
              print,'Degeneracy in positions or determinant too small.  Need more points.'
              return,map
           endif
           x3 = x[pos]
           y3 = y[pos]
           xp3 = xprime[pos]
           yp3 = yprime[pos]
           pos = pos + 1
       endif else noquit = 0

     ENDWHILE
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; find the magnitude of a
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if ( doa ) then begin
     mata = fltarr(2,2)
     mata[0,0] = r12
     mata[0,1] = r13
     mata[1,0] = yp12^2
     mata[1,1] = yp13^2

     a2 = determ( mata )/detab

     if (a2 lt 0) then print,'Correcting for negative a^2'
     a  = sqrt( abs(a2) )
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; find the magnitude of b
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if ( dob ) then begin
     matb = fltarr(2,2)
     matb[0,0] = xp12^2
     matb[0,1] = xp13^2
     matb[1,0] = r12
     matb[1,1] = r13
     
     b2 = determ( matb )/detab

     if (b2 lt 0) then print,'Correcting for negative b^2'
     b = sqrt( abs(b2) )
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;; check for reflection of axes 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if keyword_set(flipx) then a = -a
  if keyword_set(flipy) then b = -b

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;Now find costheta and sintheta
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if ( abs(detang - 0.0) lt tol ) then begin
     print,'Angle is very small'
     sintheta = 0.0
     costheta = 1.0
  endif else begin
     cosmat = fltarr(2,2)
     cosmat[0,0] = a*xp12
     cosmat[0,1] = b*yp12
     cosmat[1,0] = y12
     cosmat[1,1] = -x12
     costheta = determ( cosmat )/detang
 
     sinmat = fltarr(2,2)
     sinmat[0,0] = x12
     sinmat[0,1] = y12
     sinmat[1,0] = a*xp12
     sinmat[1,1] = b*yp12
     sintheta = determ( sinmat )/detang
  endelse

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; calculate the translations
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  c = a*xp1 - ( costheta*x1 + sintheta*y1 )
  d = b*yp1 - ( -sintheta*x1 + costheta*y1 )

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Set the map structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  map.stretch = [a,b]
  map.rot[0,0] = costheta
  map.rot[0,1] = -sintheta
  map.rot[1,0] = sintheta
  map.rot[1,1] = costheta
  map.trans = [c,d]

return,map
end











