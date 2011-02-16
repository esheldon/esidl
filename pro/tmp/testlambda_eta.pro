PRO testlambda_eta, str

  nstr = n_elements(str)
  
  d2r = !dpi/180.0d0
  r2d = 180./!dpi

  n1=200
  n2=lindgen(nstr)

  ra1 = str[n1].ra
  dec1 = str[n1].dec
  lambda1 = str[n1].lambda
  eta1 = str[n1].eta

  ra2 = str[n2].ra
  dec2 = str[n2].dec
  lambda2 = str[n2].lambda
  eta2 = str[n2].eta

  ;; old way
  t1=systime(1)

  radiff = (ra1-ra2)*d2r
  cosradiff = cos(radiff)
  sindec1 = sin(dec1*d2r)
  cosdec1 = cos(dec1*d2r)
  
  sindec2=sin(dec2*d2r)
  cosdec2=cos(dec2*d2r)

  dist = acos( sindec1*sindec2 + cosdec1*cosdec2*cosradiff )*r2d
  t1 = systime(1)-t1

  ;; new way with lambda-eta
  t2=systime(1)
  ldist = sphdist(eta1,lambda1,eta2,lambda2,/degrees)
  t2 = systime(1)-t2

;dist2 = sphdist(ra1,dec1,ra2,dec2,/degrees)
; ra -> eta in distance calculation
; dec-> lambda

  ;; old way with lambda-eta
  t3 = systime(1)
  diff = (eta1-eta2)*d2r
  cosdiff = cos(diff)
  sinlam1 = sin(lambda1*d2r)
  coslam1 = cos(lambda1*d2r)

  sinlam2 = sin(lambda2*d2r)
  coslam2 = cos(lambda2*d2r)

  ldist2 = acos( sinlam1*sinlam2 + coslam1*coslam2*cosdiff )*r2d
  t3 = systime(1)-t3

  print,t1
  print,t2
  print,t3
  print
  print, dist[30],ldist[30],ldist2[30]
  print,max( abs(dist - ldist2) )

return
END 
