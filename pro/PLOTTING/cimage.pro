PRO cimage, image, sigma, minx, maxx, miny, maxy, $
            lastpos=lastpos, resume=resume, $
            tv = tv,  $
            nx=nx, ny=ny, noneg=noneg, $
            _extra=extra

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: cimage, image [, sigma, minx, maxx, miny, maxy, lastpos=lastpos, resume=resume, tv=tv, nx=nx, ny=ny, noneg=noneg,_extra=extra]'
      return
  ENDIF 
  
  IF NOT keyword_set(tv) THEN tv = 0
  IF NOT keyword_set(noneg) THEN noneg = 0

  siz = size(image)
  sx  = siz[1]
  sy  = siz[2]

  IF n_elements(minx) NE 0 THEN BEGIN
      x = arrscl( dindgen(sx), minx, maxx )
      y = arrscl( dindgen(sy), miny, maxy )
  ENDIF 

  IF n_elements(sigma) EQ 1 THEN BEGIN 
      levels = indgen(10)+1
      IF NOT noneg THEN BEGIN 
          levels = sigma*[-reverse(levels), levels]
          c_linestyle = (levels LT 0.)
      ENDIF 
  ENDIF ELSE IF n_elements(sigma) GT 1 THEN BEGIN
      levels = sigma
      IF NOT noneg THEN BEGIN 
          levels = [-reverse(levels), levels]
          c_linestyle = (levels LT 0.)
      ENDIF 
  ENDIF ELSE BEGIN
      c_linestyle=0
      levels=0
  ENDELSE 
  c_labels = levels & c_labels[*] = 1


  IF (n_elements(nx) EQ 0) AND (n_elements(ny) EQ 0) THEN BEGIN
      nx = min([sx,sy])
      ny = nx
  ENDIF ELSE IF (n_elements(nx) EQ 0) THEN BEGIN
      nx = ny
  ENDIF ELSE IF (n_elements(ny) EQ 0) THEN BEGIN
      ny = nx
  ENDIF

  IF NOT tv THEN BEGIN 
      scrollcon, image, nx, ny, x=x, y=y, $
        levels=levels, c_linestyle=c_linestyle, c_labels=c_labels, $
        lastpos=lastpos, resume=resume, _extra=extra
  ENDIF ELSE BEGIN 
      scrolltv, image, nx, ny, _extra=extra
  ENDELSE 


return
END 
