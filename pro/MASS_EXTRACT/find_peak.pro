PRO find_peak, image, cat

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    
;       
; PURPOSE:
;    
;
; CALLING SEQUENCE:
;    
;
; INPUTS: 
;    
;
; OPTIONAL INPUTS:
;    
;
; KEYWORD PARAMETERS:
;    
;       
; OUTPUTS: 
;    
;
; OPTIONAL OUTPUTS:
;    
;
; CALLED ROUTINES:
;    
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


  IF N_params() EQ 0 THEN BEGIN 
     print,'-Syntax: find_peak, image, cat'
     print,''
     print,'Use doc_library,"find_peak"  for more help.'  
     return
  ENDIF 

;  t=systime(1)
  
  sz = size(image)
  stot = sz[4]
  sx = sz[1]
  sy = sz[2]
  
  name = 'mass1'
  ss = create_struct(name = name, $
                     'x', 0L, 'y', 0L, $
                     'val', 0., 's2n', 0., $
                     'Wild', 0L)
  
  v = replicate( image[0,0], 3, 3 )

  val = fltarr(stot)
  x = lonarr(stot)
  y = lonarr(stot)

  nobj = 0L
  FOR iy = 0L, sy-1 DO BEGIN 
      IF (iy NE 0) AND (iy NE sy-1) THEN BEGIN 
          FOR ix = 0L, sx-1 DO BEGIN 
              IF (ix NE 0) AND (ix NE sx-1) THEN BEGIN 
                  
                  v[0,0] = image[ix-1, iy-1]
                  v[1,0] = image[ix, iy-1]
                  v[2,0] = image[ix+1, iy-1]

                  v[0,1] = image[ix-1, iy]
                  v[1,1] = image[ix, iy]
                  v[2,1] = image[ix+1, iy]
                  
                  v[0,2] = image[ix-1, iy+1]
                  v[1,2] = image[ix, iy+1]
                  v[2,2] = image[ix+1, iy+1]

                  max = max(v)
                  min = min(v)

                  IF (min EQ v[1,1]) OR (max EQ v[1,1]) THEN BEGIN 
                      x[nobj] = ix
                      y[nobj] = iy
                      nobj = nobj+1
                  ENDIF 
              ENDIF 
          ENDFOR
      ENDIF
  ENDFOR 

  IF nobj EQ 0 THEN BEGIN
      print,'No objects found'
      cat = -1
  ENDIF 
  cat = replicate(ss, nobj)

  cat.x = x[0:nobj-1]
  cat.y = y[0:nobj-1]
  cat.val = image[x[0:nobj-1], y[0:nobj-1]]
                      
;  ptime,systime(1) - t

  return 
END 
