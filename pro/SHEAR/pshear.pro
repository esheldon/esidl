PRO pshear, aratio, posangle,  xsource, ysource

  IF NOT keyword_set(silent) THEN silent = 0

  ; factor to scale lenght of vector.
  factor = 1.0


  ; rescale the variables  Input must be in [0,1]

  ; make aratio start off at 1, angle in [-pi/2,pi/2] and x,y in [0,10]
  arat = 1-aratio
  theta = (2*posangle - 1)*!pi/2.
  xs = xsource*10
  ys = ysource*10
  ;; make sure r isn't too small
  xs = xs > 1.e-6
  ys = ys > 1.e-6
  r2 = xs^2 + ys^2

  xlens = 0.
  ylens = 0.

  ; find e1 and e2
  finde, arat, theta, e1, e2
  e = sqrt(e1^2 + e2^2)

  ; calculate the radial and tangential ellipticities
  e1prime = - ( e1*(xs^2 - ys^2) + e2*2.*xs*ys )/r2
  e2prime =   ( e1*2.*xs*ys  -  e2*(xs^2 - ys^2) )/r2

  ; Make the plot
  xrange=[0,10*factor]
  yrange=[0,10*factor]
  xtitle='X'
  ytitle='Y'
  ;;plot a vector from lens to source
;  plot, [xlens, xs],[ylens, ys],xrange=xrange,yrange=yrange, $
;    xtitle=xtitle,ytitle=ytitle,ystyle=1, xstyle=1
  plot, [0, 0],[0, 0],xrange=xrange,yrange=yrange, $
    xtitle=xtitle,ytitle=ytitle,ystyle=1, xstyle=1

  a = e*factor
  b = arat*a
  ; problem with angles somewhere
  tvellipse, a, b, xs, ys, theta*180/!pi, /data

  xyouts, .7*xrange[1], .95*yrange[1], 'e1 = '+ntostr(e1)
  xyouts, .7*xrange[1], .9*yrange[1], 'e2 = '+ntostr(e2)
  xyouts, .7*xrange[1], .85*yrange[1], "e1' = "+ntostr(e1prime)
  xyouts, .7*xrange[1], .8*yrange[1],  "e2' = " +ntostr(e2prime)
  xyouts, .7*xrange[1], .75*yrange[1], 'pos angle = '+ntostr(theta*180/!pi)
;  xyouts, .7*xrange[1], .7*yrange[1], 'axis ratio = '+ntostr(arat)
;  xyouts, .7*xrange[1], .65*yrange[1], 'e0 = '+ntostr(e)

  xyouts, xs+1, ys+1, '['+ntostr(xs,3)+', ' +ntostr(ys, 3)+']'

  return 
END
