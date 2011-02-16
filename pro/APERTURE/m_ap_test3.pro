PRO m_ap_test3, kappa, minx, maxx, miny, maxy, m_ap2, rsize, M_Power, lind

  IF n_params() LT 5 THEN BEGIN
      print,'-Syntax: m_ap_test3, kappa, minx, maxx, miny, maxy [, m_ap2, rsize]'
      return
  ENDIF 

  ;; Similar to m_ap_test2 but this convolves entire map with u and takes
  ;; variance.

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
  dx0 = (maxx - minx)/(sx - 1)*60. ; arcminutes

  rmax0 = 4L                    ;Pixels
  rmax_max = slow/32L

  nsize = 10              
  iii = lindgen(nsize)

  alpha = alog( float(rmax_max)/rmax0 )/(nsize - 1)
  rsize = long( rmax0*exp(alpha*iii) )
  width = 2L*rsize + 1

;print,width

  nx = sx/width
  ny = sy/width

;print
;print,nx
;print,ny


  ;; kernel sizes from 4 - slow/2.
  min_size = 2L*rmax0
  max_size = 2L*rmax_max
  
  m_ap2 = fltarr(nsize)

  sstr = ntostr(nsize)

  re = kappa
  FOR i=0, nsize-1 DO BEGIN
      print
      print,ntostr(i+1)+'/'+sstr

      u = ufunc(rsize[i])

      re = convol( temporary(re), u, /edge_truncate )

      m_ap2[i] = variance(re)

  ENDFOR 

  rsize = rsize*dx0

  w=where(m_ap2 GT 0.)
  lind = 2.16e4/rsize[w]        ;For rsize in arcminutes
  M_Power = !pi/G*m_ap2[w]/lind^2

  ptime,systime(1)-t
return
END 
