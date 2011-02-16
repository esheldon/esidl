PRO m_ap_test1, kappa, minx, maxx, miny, maxy, m_ap2, rsize, M_Power, lind

  IF n_params() LT 5 THEN BEGIN
      print,'-Syntax: m_ap_test1, kappa, minx, maxx, miny, maxy [, m_ap2, rsize]'
      return
  ENDIF 

  ;; Uses rebin to boxcar smooth the map.  Then take variance
  ;; This is alternative to m_ap_test3 which convolves map with 
  ;; smoothing function.

  L = 1.
  G = (1+L)*(2+L)^2/(1+2*L)/(3+2*L)

  ;; Sides of kappa should be powers of 2

  sz = size(kappa)
  stot = sz[4]
  sx = sz[1]
  sy = sz[2]
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; limited by smallest size
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  slow = min([sx, sy])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; slow = 2^n  n=nstep
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  nstep = long( alog(slow)/alog(2) )
  print,'Sx = ',sx
  print,'Sy = ',sy
  print,'nstep = ',nstep

  m_ap2 = fltarr(nstep)
  binsize = m_ap2

  dx0 = (maxx - minx)/(sx - 1)*60. ; arcminutes
  print,'Initial stepsize = ',dx0

  m_ap2[0] = variance(kappa)
  binsize[0] = dx0

  re = kappa
  sx = sx/2
  sy = sy/2
  dx0 = 2.*dx0
  FOR i=1, nstep-1 DO BEGIN 

      re = rebin(temporary(re), sx, sy)
      m_ap2[i] = variance(re)
      binsize[i] = dx0
      sx = sx/2
      sy = sy/2
      dx0 = 2.*dx0

  ENDFOR 

  rsize = binsize/2.

  print,'Final Stepsize = ',dx0/2.
  print

  w=where(m_ap2 GT 0.)
  lind = 2.16e4/rsize[w]        ;For rsize in arcminutes
  M_Power = !pi/G*m_ap2[w]/lind^2

return
END 
