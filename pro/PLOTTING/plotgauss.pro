PRO plotgauss, means, sigmas, $
               $
               xi_array=xi_array, $
               gauss_array=gauss_array, $
               $
               comb_xi=comb_xi, $
               comb_gauss=comb_gauss, $
               $
               overplot=overplot, $
               colors=colors, $
               _extra=_extra

  ngauss = n_elements(means)
  
  ncolor = n_elements(colors)
  IF ncolor NE 0 THEN BEGIN 
      IF ncolor EQ 1 THEN BEGIN 
          color_array = replicate(colors[0], ngauss)
      ENDIF ELSE IF ncolor EQ ngauss THEN BEGIN 
          color_array = colors
      ENDIF ELSE BEGIN 
          print,'Color must be same size as input mean/sigma or size 1'
          color_array = replicate(!p.color, ngauss)
      ENDELSE 
  ENDIF ELSE BEGIN 
      color_array = replicate(!p.color, ngauss)
  ENDELSE 


  nplot = 100
  ncomb = 1000



  xi_array = fltarr(ngauss, nplot)
  gauss_array = fltarr(ngauss, nplot)
  minx = min(means-3.5*sigmas)
  maxx = max(means+3.5*sigmas)

  comb_xi = arrscl( findgen(ncomb), minx, maxx )
  comb_gauss = fltarr(ncomb)

  FOR i=0L, ngauss-1 DO BEGIN 
      mean = means[i]
      sigma = sigmas[i]

      ;; for plotting
      txi= arrscl( findgen(nplot), mean - 3.5*sigma, mean+3.5*sigma )
      tgauss = gaussprob(txi, mean, sigma)

      xi_array[i, *] = txi
      gauss_array[i, *] = tgauss

      ;; for combined
      tgauss = gaussprob(comb_xi, mean, sigma)
      comb_gauss[*] = comb_gauss[*] + tgauss

  ENDFOR 

  ;; normalize combined
  comb_gauss = comb_gauss/qgauss(comb_gauss, comb_xi, 1000)

  maxy = max(gauss_array)
  
  IF n_elements(xrange) EQ 0 THEN xrange = [minx, maxx]
  IF n_elements(yrange) EQ 0 THEN yrange = [0, maxy]

  FOR i=0L, ngauss-1 DO BEGIN 

      color = color_array[i]

      xi = xi_array[i,*]
      gauss = gauss_array[i, *]
      IF i EQ 0 THEN BEGIN 
          IF keyword_set(overplot) THEN BEGIN 
              oplot, xi, gauss,color=color, _extra=_extra
          ENDIF ELSE BEGIN 
              plot, xi, gauss, xrange=xrange, yrange=yrange,color=color, $
                _extra=_extra
          ENDELSE 
      ENDIF ELSE BEGIN 
          oplot, xi, gauss,color=color, _extra=_extra
      ENDELSE 
  ENDFOR 

  
  


END 
