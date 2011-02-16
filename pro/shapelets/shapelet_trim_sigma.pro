FUNCTION shapelet_trim_sigma, image, center, sigma, nsig, minx=minx, maxx=maxx, miny=miny, maxy=maxy

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: new = shapelet_trim_sigma(image, center, sigma, nsig, minx=, maxx=, miny=, maxx=)'
      return,-1
  ENDIF 

  ;; trim 3.5 to sigma region
  imsize = size(image, /dim)

  minx = (center[0] - nsig[0]*sigma[0]) > 0
  maxx = (center[0] + nsig[0]*sigma[0]) < (imsize[0]-1)
  miny = (center[1] - nsig[0]*sigma[0]) > 0
  maxy = (center[1] + nsig[0]*sigma[0]) < (imsize[1]-1)
  
  IF ((minx GT 0) AND (maxx LT (imsize[0]-1)) AND $
      (miny GT 0) AND (maxy LT (imsize[1]-1)) ) THEN BEGIN 
      
      new_image = image[minx:maxx, miny:maxy]
      return,new_image      
  ENDIF ELSE BEGIN 
      return,image
  ENDELSE 

END 
