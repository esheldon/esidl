PRO test_circ_rect_intersect

  !p.multi=[0,2,4]
  pcold=!p.charsize
  !p.charsize=0.75

  ;; set up region

  x1=0.0
  y1=0.0

  x2=1.0
  y2=0.0

  x3=1.0
  y3=1.0

  x4=0.0
  y4=1.0
  
  ;; some points for monte carlo
  ntot = 1.e6
  x=arrscl( randomu(seed, ntot), x1, x2, arrmin=0., arrmax=1. )
  y=arrscl( randomu(seed, ntot), y1, y3, arrmin=0., arrmax=1. )
  areatot=(x2-x1)*(y3-y1)
  density = ntot/areatot

  ;;;;;;;;;;;;;;;;
  ;; no cross
  ;;;;;;;;;;;;;;;;

  print,'No crossing'
  cenx=0.5
  ceny=0.5
  radius=0.25

  dfromright = x2-cenx
  dfromleft = cenx-x1
  dfromtop = y3-ceny
  dfrombottom = ceny-y1

  area = !pi*radius^2
  print,'truth: ',area
  print,'----------------------------------------'
  circ_rect_intersect,dfromtop, dfrombottom, dfromright, dfromleft,radius,carea
  print,carea

  aplot,1,[0],xrange=[-0.1, 1.1], yrange=[-0.1,1.1],$
    title='Area = '+ntostr(area)+'  Measured = '+ntostr(carea)
  plot_box, x1, x2, y1, y3
  tvcircle, radius, cenx, ceny,/data

  if display_type() eq 'X' then key=get_kbrd(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; entire rectangle encompassed
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,'Entire rectanble encompassed'
  cenx=0.5
  ceny=0.5
  radius=0.6*sqrt(2.)

  dfromright = x2-cenx
  dfromleft = cenx-x1
  dfromtop = y3-ceny
  dfrombottom = ceny-y1

  area = areatot
  print,'truth: ',area
  print,'----------------------------------------'
  circ_rect_intersect,dfromtop, dfrombottom, dfromright, dfromleft,radius,carea
  print,carea

  aplot,1,[0],xrange=[-0.1, 1.1], yrange=[-0.1,1.1],$
    title='Area = '+ntostr(area)+'  Measured = '+ntostr(carea)
  plot_box, x1, x2, y1, y3
  tvcircle, radius, cenx, ceny,/data

  if display_type() eq 'X' then key=get_kbrd(1)


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; cross top/bottom/right/left
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  print,'cross one side'

  radius=0.25

  nvals=4
  xvals = [0.5, 0.5, 0.8, 0.2]
  yvals = [0.2, 0.8, 0.5, 0.5]

  area = binarea(y3-yvals[1], 0.0, radius)
  print,'truth: ',area
  print,'----------------------------------------'

  
  FOR i=0L, nvals-1 DO BEGIN 
      cenx=xvals[i]
      ceny=yvals[i]
      dfromright = x2-cenx
      dfromleft = cenx-x1
      dfromtop = y3-ceny
      dfrombottom = ceny-y1

      circ_rect_intersect,dfromtop, dfrombottom, dfromright, dfromleft,radius,carea
      print,carea

      aplot,1,[0],xrange=[-0.1, 1.1], yrange=[-0.1,1.1],$
        title='Area = '+ntostr(area)+'  Measured = '+ntostr(carea)
      plot_box, x1, x2, y1, y3
      tvcircle, radius, cenx, ceny,/data

      if display_type() eq 'X' then key=get_kbrd(1)
  ENDFOR 



  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; corner outside
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; must use monte carlo to test
  print,'Corner outside circle'

  radius=0.25
 
  nval=4
  xvals=[0.1,0.9,0.9,0.1]
  yvals=[0.1,0.1,0.9,0.9]

  ;; same for all, just do one
  w=where(sqrt((xvals[0]-x)^2 + (yvals[0]-y)^2) LT radius, nw)
  area = nw/density

  print,'truth: ',area
  print,'----------------------------------------'

  FOR i=0L, nval-1 DO BEGIN 
      cenx=xvals[i]
      ceny=yvals[i]
      
      dfromright = x2-cenx
      dfromleft = cenx-x1
      dfromtop = y3-ceny
      dfrombottom = ceny-y1

      circ_rect_intersect,dfromtop, dfrombottom, dfromright, dfromleft,radius,carea
      print,carea

      aplot,1,[0],xrange=[-0.1, 1.1], yrange=[-0.1,1.1],$
        title='Area (MC) = '+ntostr(area)+'  Measured = '+ntostr(carea)
      plot_box, x1, x2, y1, y3
      tvcircle, radius, cenx, ceny,/data

      if display_type() eq 'X' then key=get_kbrd(1)
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;
  ;; corner inside
  ;;;;;;;;;;;;;;;;;;;;;

  ;; must use monte carlo to test
  print,'Corner outside circle'

  radius=0.25
 
  nval=4
  xvals=[0.2,0.8,0.8,0.2]
  yvals=[0.2,0.2,0.8,0.8]

  ;; same for all, just do one
  w=where(sqrt((xvals[0]-x)^2 + (yvals[0]-y)^2) LT radius, nw)
  area = nw/density

  print,'truth: ',area
  print,'----------------------------------------'

  FOR i=0L, nval-1 DO BEGIN 
      cenx=xvals[i]
      ceny=yvals[i]
      
      dfromright = x2-cenx
      dfromleft = cenx-x1
      dfromtop = y3-ceny
      dfrombottom = ceny-y1

      circ_rect_intersect,dfromtop, dfrombottom, dfromright, dfromleft,radius,carea
      print,carea

      aplot,1,[0],xrange=[-0.1, 1.1], yrange=[-0.1,1.1],$
        title='Area (MC) = '+ntostr(area)+'  Measured = '+ntostr(carea)
      plot_box, x1, x2, y1, y3
      tvcircle, radius, cenx, ceny,/data

      if display_type() eq 'X' then key=get_kbrd(1)
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;
  ;; All corners outside
  ;;;;;;;;;;;;;;;;;;;;;

  ;; must use monte carlo to test
  print,'All corners outside

  radius=0.6
 
  nval=4
  cenx=0.5
  ceny=0.5

  ;; same for all, just do one
  w=where(sqrt((cenx-x)^2 + (ceny-y)^2) LT radius, nw)
  area = nw/density

  print,'truth: ',area
  print,'----------------------------------------'
      
  dfromright = x2-cenx
  dfromleft = cenx-x1
  dfromtop = y3-ceny
  dfrombottom = ceny-y1

  circ_rect_intersect,dfromtop, dfrombottom, dfromright, dfromleft,radius,carea
  print,carea

  aplot,1,[0],xrange=[-0.1, 1.1], yrange=[-0.1,1.1],$
    title='Area (MC) = '+ntostr(area)+'  Measured = '+ntostr(carea)
  plot_box, x1, x2, y1, y3
  tvcircle, radius, cenx, ceny,/data

  !p.multi=0
  !p.charsize=pcold

END 
