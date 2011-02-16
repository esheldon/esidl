PRO voronoi_simulate, cat, dens, simcat, simdens


  n=n_elements(cat)
  IF n_elements(simcat) EQ 0 THEN BEGIN 
      maxx = max(cat.ra)
      minx = min(cat.ra)
  
      maxy = max(cat.dec)
      miny = min(cat.dec)


      sscat = create_struct('RA', 0D, 'DEC', 0D)

      simcat = replicate(sscat, n)

      simcat.ra = arrscl( randomu(seed,n), minx, maxx, arrmin=0., arrmax=1.)
      simcat.dec = arrscl( randomu(seed,n), miny, maxy, arrmin=0., arrmax=1.)

      voronoi_density, cat.ra, cat.dec, dens
      voronoi_density, simcat.ra, simcat.dec, simdens
  ENDIF 

;  pold = !p.multi
;  !p.multi = [0,1,2]
  
  t1 = 'Run 752/756 foreground gals'
  t2 = 'Random'
  position = [1.e4,800]

  tit=ntostr(n)+' gals/simulated points'
  yt = 'N'
  xt = 'Density (arb. units)'
  bin=500
  yrange=[0,1000.]
  xrange=[-500,2.e4]


  plothist, dens, bin=bin, yrange=yrange, xrange=xrange, $
            xtit=xt,ytit=yt,tit=tit
  plothist,simdens,bin=bin,yrange=yrange,xrange=xrange, /overplot, linestyle=2

  legend,[t1,t2],linestyle=[0,2],position=position

;  !p.multi=pold
  return
END 
