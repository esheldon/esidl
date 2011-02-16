PRO m_ap_test2, kappa, minx, maxx, miny, maxy, m_ap2, rsize, M_Power, lind, m_ap2err, M_Powererr

  IF n_params() LT 5 THEN BEGIN
      print,'-Syntax: m_ap_test2, kappa, minx, maxx, miny, maxy [, m_ap2, rsize]'
      return
  ENDIF 

  ;; Finds a set of "independent" points in the map and convolves with
  ;; smoothing function u to get m_ap at those points.  Then take variance.
  ;; Repeat for vaying smoothing lengths.

  t=systime(1)

  L = 1.
  G = (1+L)*(2+L)^2/(1+2*L)/(3+2*L)

  sz = size(kappa)
  stot = sz[4]
  sx = sz[1]
  sy = sz[2]
  
  index = lindgen(stot)
  xind = index/sx
  yind = index MOD sx


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; limited by smallest size
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  slow = min([sx, sy])
  dx0 = (maxx - minx)/(sx - 1)*60. ; pixel size in arcminutes

  rmax0 = 4L                    ;Pixels
  rmax_max = slow/8L

  nsize = 10              
  iii = lindgen(nsize)

  alpha = alog( float(rmax_max)/rmax0 )/(nsize - 1)
  rsize = long( rmax0*exp(alpha*iii) )
  width = 2L*rsize + 1

print,width

  nx = sx/width
  ny = sy/width

  ;; kernel sizes from 4 - slow/2.
  min_size = 2L*rmax0
  max_size = 2L*rmax_max
  
  m_ap2 = fltarr(nsize)
  m_ap2err = m_ap2

  sstr = ntostr(nsize)

  
  FOR i=0, nsize-1 DO BEGIN
      print
      print,ntostr(i+1)+'/'+sstr,'  ndots= ',ntostr(ny[i])

      val = fltarr( nx[i]*ny[i] )
      u = ufunc(rsize[i])

      count = 0L

      FOR iy=0L, ny[i]-1 DO BEGIN
          print,format='($,A)','.'

          y = rsize[i] + iy*width[i]
          ydiff = abs(y - yind)

          w1 = where( ydiff LE rsize[i], nw1)
          FOR ix=0L, nx[i]-1 DO BEGIN
              x = rsize[i] + ix*width[i]
              xdiff = abs(x - xind[w1])
              
              w2 = where( xdiff LE rsize[i], nw2)
              w=w1[w2]

              val[count] = total( u*kappa[w] )
              
              count = count + 1
          ENDFOR 
      ENDFOR 
      print
      m_ap2[i] = variance(val)
      m_ap2err[i] = m_ap2[i]/sqrt(2.*nx[i]*ny[i])

      print,ntostr(m_ap2[i]),' +/- ',ntostr(m_ap2err[i])

  ENDFOR 

  rsize = rsize*dx0

  w=where(m_ap2 GT 0.)
  lind = 2.16e4/rsize[w]        ;For rsize in arcminutes
  M_Power = !pi/G*m_ap2[w]/lind^2
  M_Powererr = !pi/G*m_ap2err[w]/lind^2

  ptime,systime(1)-t
return
END 
