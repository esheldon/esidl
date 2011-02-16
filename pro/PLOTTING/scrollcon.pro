PRO scrollcon, image, nx, ny, x=x, y=y, $
               lastpos=lastpos, resume=resume, _extra=extra

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: scrollcon, image [, nx, ny, x=x, y=y, lastpos=lastpos, resume=resume, _extra=extra]'
      return
  ENDIF 
  IF NOT keyword_set(resume) THEN resume=0
  IF (n_elements(x) NE 0) AND (n_elements(y) NE 0) THEN doaxis=1 ELSE doaxis=2
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

  IF resume AND (n_elements(lastpos) NE 0) THEN BEGIN 
      ix = lastpos[0]
      iy = lastpos[1]
  ENDIF ELSE BEGIN 
      ix=0
      iy=0
  ENDELSE 
  ixold=ix
  iyold=iy
  WHILE 1 DO BEGIN
      x2 = ix + bx-1
      y2 = iy + by-1
      IF (x2 LE sx-1) AND (ix NE -1) AND $
         (y2 LE sy-1) AND (iy NE -1) THEN BEGIN

          IF doaxis THEN BEGIN 
              contour, image[ix:x2,iy:y2], x[ix:x2], y[iy:y2], $
                       xstyle=1, ystyle=1, _extra=extra

          ENDIF ELSE BEGIN 
              contour, image[ix:x2,iy:y2], xstyle=1, ystyle=1, $
                       _extra=extra
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
              '': BEGIN & lastpos=[ix,iy] & return & END 
              ELSE: go=1
          ENDCASE 
      ENDWHILE 
  ENDWHILE 

  return
END 
          
