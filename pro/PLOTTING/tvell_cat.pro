pro tvell_cat, cat, cattype, cindex=cindex, scalef=scalef,color=color,$
            range=range, step=step, adaptive=adaptive

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: 
;    TVELL_CAT
;       
; PURPOSE: 
;    wrapper for tvellipse.  Takes a Photo or SExtractor catalog file
;    and, for each object, draws an ellipse on current plot window.
;	   
;
; CALLING SEQUENCE: tvell_cat, cat, cattype, cindex=cindex, scalef=scalef,
;           color=color,range=range
;      
;                 
;
; INPUTS:  cat: Photo or SExtractor catalog.
;          cattype: cattype = 1 for photo catalog.
;                   cattype = 2 for SExtractor catalog.
;
; OPTIONAL INPUTS:
;           cindex:  The color index to use from the photo catalog.  This 
;                    assumes photo convention 0..4 -> u,g,r,i,z
;           scalef:  Factor by which to multiply the ellipses axes.
;           color:   Color of ellipse.  0-!d.n_colors  (black-white)
;           range:   Range of plot in which the object must be to get
;                    an ellipse.
;
; KEYWORD PARAMETERS:
;           /step: step though and circle one at a time with prompting.
;           /adaptive: if set, it assumes there are adaptive moment
;                      e1 and e2 in sextractor catalog named e1_ad, e2_ad.
;                      same with x2_ad, y2_ad for moments.
;
; CALLED ROUTINES: TVELLIPSE
; 
; PROCEDURE: Calculate ellipse parameters based on size and shape information
;            in each catalog.  Find position of each object.  Call
;            TVELLIPSE for each object.
;	
;	
;
; REVISION HISTORY: Author  Erin Scott Sheldon  UM  3/23/99
;	
;       
;                                      
;-                                        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF N_params() eq 0 then begin
     print,'-Syntax: tvell_cat, cat, cattype, scalef=scalef, color=color,range=range, step=step'
     print,''
     print,'cattyp = 1 for photo,  2 for SExtractor'
     return
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; some parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(scalef) EQ 0 THEN scalef = 1.0
  IF n_elements(cindex) EQ 0 THEN cindex=2
  IF n_elements(color) EQ 0 THEN color = 0.9*!d.n_colors
  IF n_elements(range) EQ 0 THEN range=[[0.0,1.0e20],[0.0,1.0e20]]
  IF keyword_set(step) THEN step = 1 ELSE step = 0
  IF keyword_set(adaptive) THEN adaptive = 1 ELSE adaptive = 0

  nn = n_elements(cat) 
  maxsize = 500.                ;pixels
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Set up the tvellipse inputs for the input catalog
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                ; setup for photo
  IF cattype EQ 1 THEN BEGIN 
    x=cat.colc[cindex]
    y=cat.rowc[cindex]
    a = sqrt( cat.iso_a[cindex]^2 + cat.iso_b[cindex]^2)

    qu2e, cat.q[cindex], cat.u[cindex],e1,e2
    
    w=where( sqrt(e1^2 + e2^2) LE 1.0 AND $
             abs(a) LT maxsize AND $
             x GE range[0,0] AND x le range[1,0] AND $
             y GE range[0,1] AND y LE range[1,1] , nw)

  ENDIF ELSE IF cattype EQ 2 THEN BEGIN ; setup SExtractor catalog

    x=cat.x_image
    y=cat.y_image
    IF adaptive THEN BEGIN 
      a = sqrt( cat.x2_ad^2+ cat.y2_ad^2 )
      e1 = cat.e1_ad
      e2 = cat.e2_ad
    ENDIF ELSE BEGIN 
      a = sqrt( cat.x2_image^2+ cat.y2_image^2 )
      e1 = cat.e1
      e2 = cat.e2
    ENDELSE 

    w=where( sqrt(e1^2 + e2^2) LE 1.0 AND $
             abs(a) LT maxsize AND $
             x GE range[0,0] AND x LE range[1,0] AND $
             y GE range[0,1] AND y LE range[1,1], nw)

  ENDIF ELSE BEGIN
    print,'Invalid Catalog Type'
    return
  ENDELSE 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Draw an ellipse on each object
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF nw NE 0 THEN BEGIN 
    a=scalef*a
    FOR jj = 0, nw-1 DO BEGIN 
      i = w[jj]
      findabtheta, e1[i], e2[i], aratio, pos_ang,/lupton
      pos_ang = pos_ang*180./!pi
      b = aratio*a[i]
         
      tvellipse, a[i], b, x[i], y[i], pos_ang, color, /data
      IF (keyword_set(step)) THEN BEGIN
        print,'Object index: ',ntostr(i)
        print,'Hit a key for next object   (q to quit)'
        key = get_kbrd(20)
        IF (key EQ 'q') THEN return
      ENDIF
    ENDFOR 
  ENDIF 



  return 
END 






