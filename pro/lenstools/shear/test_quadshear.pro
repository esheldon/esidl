
PRO test_quadshear, struct, rotate=rotate, g=g

  radius = 4.0                  ;arcmin
  radius = radius*60.           ;arcsec
  radius = radius/0.4           ;pixels

  run=756
  camcol=3
  rerun=1
  field=300

  d2r=!dpi/180.
  r2d=1d/d2r

  fetch_dir,run,camcol,rerun,dir,atldir,corrdir=corrdir,corratldir=corratldir

  IF n_elements(struct) EQ 0 THEN BEGIN
      
      read_tsobj,corrdir,st,start=field-1,nf=3,front='adatc'
      add_tag,st,'fibercounts',fltarr(5),struct

      struct.fibercounts = struct.petrocounts

  ENDIF 

  w=where((struct.field EQ field) AND $
          (struct.e1[2] NE 1.e10) AND $
          (struct.petrocounts[2] LT 18) AND $
          (struct.petrocounts[2] GT 16) AND $
          (struct.type[1] EQ 3) AND $
          (struct.type[2] EQ 3) )

  help,w

  ;;;;;;;;;;;;;;;;
  ;; ra/dec
  ;;;;;;;;;;;;;;;

  cendec = struct[w[g]].dec
  cenra  = struct[w[g]].ra

  radiff = (cenra-struct.ra)*d2r
  cosradiff = cos(radiff)
                  
  tcendec = cendec*d2r
  sincendec=sin(tcendec)
  coscendec=cos(tcendec)

  sinscatdec = sin(struct.dec*d2r)
  cosscatdec = cos(struct.dec*d2r)

  ;; distance in degrees
  Rradec=acos( sincendec*sinscatdec + $
               coscendec*cosscatdec*cosradiff )*r2d*3600.

  thetaradec=atan(sin(radiff), $
                  (sincendec*cosradiff - $
                   coscendec*sinscatdec/cosscatdec) )-!dpi
  wwradec=where(Rradec LE radius*0.4)
  xrelradec = Rradec[wwradec]*cos(thetaradec)
  yrelradec = Rradec[wwradec]*sin(thetaradec)

  ;;;;;;;;;;;;;;;
  ;; lambda/eta
  ;;;;;;;;;;;;;;;
  ;; copying from objshear_lambda
  eq2survey, struct.ra,struct.dec,lambda,eta
  ceneta = eta[w[g]]
  cenlam = lambda[w[g]]
print,'Ceneta = ',ceneta,'  Cenlam = ',cenlam

  etadiff = (ceneta-eta)*d2r
  cosetadiff = cos(etadiff)

  tcenlam = cenlam*d2r          ;lens center in radians
  sincenlam = sin(tcenlam)
  coscenlam = cos(tcenlam)

  tscatlam = lambda*d2r
  sinscatlam = sin(tscatlam)
  cosscatlam = cos(tscatlam)

  args = sincenlam*sinscatlam + coscenlam*cosscatlam*cosetadiff
  warg = where(args LT -1., narg)
  IF narg NE 0 THEN args[warg] = -1.
  warg = where(args GT 1., narg)
  IF narg NE 0 THEN args[warg] = 1.
                  
  R=acos( args )*r2d*3600.

  ww=where(R LE radius*0.4)

  theta=atan(sin(etadiff[ww]), $
             (sincenlam*cosetadiff[ww] - $
              coscenlam*sinscatlam[ww]/cosscatlam[ww]) ) -!dpi

  xrel=  R[ww]*sin(theta)
  yrel=  R[ww]*cos(theta)

  e1 = struct[w[g]].e1[2]
  e2 = struct[w[g]].e2[2]
  IF keyword_set(rotate) THEN BEGIN
      rotate_e1e2, struct[w[g]].rotation[2], e1, e2,$
        e1out, e2out
      e1=e1out
      e2=e2out
  ENDIF 

  e1old=e1
  e2old=e2
  ltheta=0.5*atan(e2,e1)

  print,'Angle = ',ltheta*r2d
  rotate_xy, ltheta, xrel, yrel, newx, newy

  !p.multi=[0,0,2]
  xmm=median(xrel)
  xrange=[xmm-250,xmm+250]
  yrange=[-250,250]
  aplot,1,xrel,yrel,psym=4,xrange=xrange,yrange=yrange
;aplot,1,xrelradec,yrelradec,psym=4,yrange=yrange
  aplot,1,(struct[ww].dec-cendec)*3600., (struct[ww].ra-cenra)*3600.,psym=4
  
  key=get_kbrd(1)
  IF key EQ 'q' THEN return

  !p.multi=[0,2,2]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; First all neighbors, then in half space
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  aplot,1,xrel,yrel,psym=3,xrange=xrange,yrange=yrange
  tvellipse,20.,10.,0.0,0.0,ltheta*r2d,/data

  IF keyword_set(rotate) THEN BEGIN 
      rotate_e1e2, struct[ww].rotation[2], struct[ww].e1[2], struct[ww].e2[2],$
        e1neigh, e2neigh

  ENDIF ELSE BEGIN 
      e1neigh = struct[ww].e1[2]
      e2neigh = struct[ww].e2[2]
  ENDELSE 

  angle=0.5*atan(e2neigh,e1neigh)
  nobj = n_elements(angle)
  FOR i=0L, nobj-1 DO tvellipse,10, 5, xrel[i], yrel[i], $
    angle[i]*r2d, /data

  rotate_e1e2, ltheta, e1neigh, e2neigh, $
	ne1, ne2
  rotate_e1e2, ltheta, e1, e2, newle1, newle2
  e1=newle1
  e2=newle2
  ltheta=0.5*atan(e2,e1)

  angle=0.5*atan(ne2,ne1)

  aplot, 1, newx, newy, psym=3, xrange=xrange,yrange=yrange
  tvellipse,20.,10.,0.0,0.0,ltheta*r2d,/data
  FOR i=0L, nobj-1 DO tvellipse, 10, 5, newx[i], newy[i], $
    angle[i]*r2d, /data


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; In a half space
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  half1 = where( abs(newy) GT abs(newx), num1) 
  half2 = where( abs(newx) GE abs(newy), num2) ;along major axis

  half=half1
  ww = ww[half]

  xrel=xrel[half] & yrel=yrel[half]
  newx=newx[half] & newy=newy[half]
  e2neigh=e2neigh[half] & e1neigh=e1neigh[half]
  angle=0.5*atan(e2neigh,e1neigh)
  nobj=n_elements(half)

  e1=e1old
  e2=e2old
  ltheta=0.5*atan(e2,e1)
  aplot,1,xrel,yrel,psym=3,xrange=xrange,yrange=yrange
  tvellipse,20.,10.,0.0,0.0,ltheta*r2d,/data

  FOR i=0L, nobj-1 DO tvellipse,10, 5, xrel[i], yrel[i], $
    angle[i]*r2d, /data

  rotate_e1e2, ltheta, e1neigh, e2neigh, $
	ne1, ne2
  rotate_e1e2, ltheta, e1, e2, newle1, newle2
  e1=newle1
  e2=newle2
  ltheta=0.5*atan(e2,e1)

  angle=0.5*atan(ne2,ne1)

  aplot, 1, newx, newy, psym=3, xrange=xrange,yrange=yrange
  oplot,[-1000,1000],[-1000,1000]
  oplot,[-1000,1000],[1000,-1000]
  tvellipse,20.,10.,0.0,0.0,ltheta*r2d,/data
  FOR i=0L, nobj-1 DO tvellipse, 10, 5, newx[i], newy[i], $
    angle[i]*r2d, /data

  !p.multi=0

  return

  key=get_kbrd(1)
  IF key EQ 'q' THEN return

;  get_atlas,struct,w,dir=atldir

  ra = [struct[w[g]].ra, struct[ww].ra ]
  dec = [struct[w[g]].dec, struct[ww].dec]

  fchart_circ_radec,struct,ra,dec,fchart,$
    clr=2,radius=radius,circ_rad=20

  print,'Angle = ',ltheta*r2d

END 
