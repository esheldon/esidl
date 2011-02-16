FUNCTION shapelet_trim_image, image, sky, skysig, nskysig, $
                              smoothing_window=smoothing_window, $
                              slack=slack, $
                              minx=minx, maxx=maxx, $
                              miny=miny, maxy=maxy, $
                              smoothed_image=smoothed_image, $
                              status=status

  on_error, 2
  status = 1

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: new = shapelet_trim_image(image, sky, skysig, nskysig, smoothing_window=, slack=, minx=, maxx=, miny=, maxy=, smoothed_image=,status=)'
      return,-1
  ENDIF 

  ;; The smoothing window
  IF n_elements(window) EQ 0 THEN window=4

  ;; Get image dimensions
  sz=size(image)

  nx = sz[1]
  ny = sz[2]
  ntot = sz[4]

  ;; A 2-d index into the image
  index = lindgen(ntot)
  x = index MOD nx
  y = index/nx

  ;; Smooth the image. /edge_truncate is crucial
  timage = smooth(image, window, /edge_truncate)

  w=where(timage GE sky+nskysig*skysig, nw)

  IF nw EQ 0 THEN return, -1


  IF arg_present(smoothed_image) THEN BEGIN 
      smoothed_image=temporary(timage)
  ENDIF ELSE BEGIN 
      timage = 0
  ENDELSE 

  ;; Get the boundary
  minx = min(x[w], max=maxx)
  miny = min(y[w], max=maxy)

  ;; Add slack pixels
  IF n_elements(slack) NE 0 THEN BEGIN 
      minx = minx - slack[0] > 0
      maxx = maxx + slack[0] < nx-1

      miny = miny - slack[0] > 0
      maxy = maxy + slack[0] < ny-1
  ENDIF 

  new_image=image[minx:maxx, miny:maxy]

  status = 0
  return, new_image

END 
