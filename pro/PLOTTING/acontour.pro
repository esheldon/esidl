PRO acontour, aspect, z, x, y, center=center, _extra=extra, xtick_get=xtick_get, ytick_get=ytick_get

  np = n_params()
  IF (np LT 2) OR (np GT 4) THEN BEGIN 
     print,'-Syntax: acontour, aspect, z, [x, y,] center=center,_extra=extra'
     print,''
     print,' aspect = xsize/ysize'
     print,'Use doc_library,"aplot"  for more help.'  
     return
  ENDIF 

  IF !p.multi[1] EQ 0 THEN !p.multi[1] = 1
  IF !p.multi[2] EQ 0 THEN !p.multi[2] = 1

  plot, [0,1], [0,1], /nodata, xstyle=4, ystyle=4
  px = !x.window*!d.x_vsize
  py = !y.window*!d.y_vsize

  xsize = px[1] - px[0]
  ysize = py[1] - py[0]

  CASE 1 OF 
      (aspect EQ 1): IF xsize GT ysize THEN xsize=ysize ELSE ysize = xsize  
      (aspect LT 1): xsize = ysize*aspect
      (aspect GT 1): ysize = xsize/aspect
  ENDCASE 

  px[1] = px[0] + xsize
  py[1] = py[0] + ysize

  IF keyword_set(center) THEN dcenter, xsize, ysize, px, py

  position = [ [ px(0), py(0)], [ px(1), py(1) ] ]

  IF np EQ 4 THEN BEGIN         ;Case where x is really y
      contour, z, x, y, position=position, /device, /noerase, _extra=extra, xtick_get=xtick_get, ytick_get=ytick_get
  ENDIF ELSE IF np EQ 3 THEN BEGIN 
      contour, z, x, position=position, /device, /noerase, _extra=extra, xtick_get=xtick_get, ytick_get=ytick_get
  ENDIF ELSE BEGIN
      contour, z, position=position, /device, /noerase, _extra=extra, xtick_get=xtick_get, ytick_get=ytick_get
  ENDELSE 

END 
