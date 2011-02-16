PRO scrolltv, image, nx, ny, x=x, y=y, _extra=extra

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: scrolltv, image [, nx, ny, x=x, y=y, _extra=extra]'
      return
  ENDIF 
  IF n_elements(x) NE 0 THEN doaxis=1 ELSE doaxis=2
  siz = size(image)
  sx = siz[1]
  sy = siz[2]

  IF n_elements(nx) EQ 0 THEN bx = sx ELSE bx = nx < sx
  IF n_elements(ny) EQ 0 THEN by = sy ELSE by = ny < sy

  
  up   = '[A'
  down = '[B'
  left = '[D'
  right = '[C'
  pgup='[5~'
  pgdown = '[6~'

  ix=0
  iy=0
  ixold=ix
  iyold=iy
  WHILE 1 DO BEGIN
      x2 = ix + bx-1
      y2 = iy + by-1
      IF (x2 LE sx-1) AND (ix NE -1) AND $
         (y2 LE sy-1) AND (iy NE -1) THEN BEGIN

          IF doaxis THEN BEGIN 
              mrdis, image[ix:x2,iy:y2], _extra=extra, nsig=2.0, /silent

          ENDIF ELSE BEGIN 
              mrdis, image[ix:x2,iy:y2], _extra=extra, nsig=2.0, /silent
          ENDELSE 

      ENDIF ELSE BEGIN
          ix = ixold
          iy = iyold
      ENDELSE 
      go=1
      WHILE go DO BEGIN
          getseq, key
          CASE key OF 
              up: BEGIN & iyold=iy & iy=iy+1 & go=0 & END
              down: BEGIN & iyold=iy & iy=iy-1 & go=0 & END 
              right: BEGIN & ixold=ix & ix=ix+1 & go=0 & END
              left:  BEGIN & ixold=ix & ix=ix-1 & go=0 & END
              pgup:  BEGIN & iyold=iy & iy=iy+by-1 < (sy-by) & go=0 & END 
              pgdown: BEGIN & iyold=iy & iy=iy-(by-1) > 0 & go=0 & END 
              '': return
              ELSE: go=1
          ENDCASE 
      ENDWHILE 
  ENDWHILE 

  return
END 
          
