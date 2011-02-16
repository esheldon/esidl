PRO implot_setup, image, xsize, ysize, px, py, xrng, yrng, center=center

  ;; Here is how I use this:

  sz = size(image, /dimensions)

  nx = sz[0]
  ny = sz[1]

  aspect = float(nx)/ny

  plot, [0,1],[0,1],/nodata,xstyle=4,ystyle=4
  px=!x.window*!d.x_vsize
  py=!y.window*!d.y_vsize
  xsize=px[1]-px[0]
  ysize=py[1]-py[0]
  
  IF xsize GT ysize*aspect THEN xsize=ysize*aspect ELSE ysize=xsize/aspect 
  px[1]=px[0]+xsize
  py[1]=py[0]+ysize

  nxm=nx-1
  nym=ny-1

  IF n_elements(xrange) EQ 0 THEN BEGIN
      xrng=[ -0.5, nxm+0.5]
  ENDIF ELSE BEGIN
      xrng=[xrange(0), xrange(n_elements(xrange)-1)]
  ENDELSE 

  IF n_elements(yrange) EQ 0 THEN BEGIN
      yrng = [-0.5,nym+0.5]
  ENDIF ELSE BEGIN
      yrng = [yrange(0), yrange(n_elements(yrange)-1)]
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; center up the display
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(center) THEN dcenter, xsize, ysize, px, py, /silent


END 
