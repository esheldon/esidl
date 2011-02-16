PRO myplot, vdisp, Zlens, Zsource

  IF n_params() EQ 0 THEN BEGIN 
    print,'-Syntax:  myplot, ....'
    return
  ENDIF

  ;; rescale variables
  vdscale = 440.0               ;km/s
  vd = vdisp*vdscale
  Zl = (Zlens > .00001)/2.0
  Zs = (Zsource > .00001)/2.0

  n=100
  minr = 25.0
  maxr=1000.0

  radius = arrscl( findgen(n), minr, maxr)

  xrange=[0,max(radius)]
  xtitle='Radius (kpc)'
  ytitle='Shear'

  shear = shearsis(vd, Zs, Zl, kpc=radius)
  help,shear
print,max(shear)
  yrange=[0,.15]
  plot, radius, shear,xtitle=xtitle,ytitle=ytitle,xrange=xrange, $
    yrange=yrange,xstyle=1,ystyle=1

  xyouts, 400, .95*yrange[1],'vdisp = '+ntostr(vd)+' km/s'
  xyouts, 400, .85*yrange[1],'Z lens = '+ntostr(Zl)
  xyouts, 400, .75*yrange[1],'Z source = '+ntostr(Zs)



return
END


