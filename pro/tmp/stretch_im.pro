PRO stretch_im


aplot, 1.0, [-1,1],[-1,1],yrange=[-1,1],xrange=[-1,1], /nodata, ystyle=4, xstyle=4

e1=-.005
e2= 0.
            
findabtheta, e1, e2, aratio, theta
            
area = .01
a = sqrt(area/aratio)
b=aratio*a

;xi=.23
;yi=0.
;tvellipse, a, b, xi,yi, theta*180./!pi, /data
 
rad=sqrt(area/!pi)*1.7
xs=0.
ys=0.
tvcircle, rad, xs,ys, /data

xlens = -.6
ylens = 0.
tvcircle, .1*rad, xlens,ylens, /data

xi=0.01
yi=0.
tvellipse, a, b, xi,yi, theta*180./!pi, /data

fac=1.55
;arrow, -.02, .98*rad, xi, fac*rad,/data,Hsize = !D.X_SIZE /80.
;arrow, -.02, -.98*rad, xi, -fac*rad,/data,Hsize = !D.X_SIZE /80.


xyouts, 1.2*xlens, -.55*rad, 'Lens', /data
xyouts, -0.1, -2.*rad, 'Source', /data
xyouts, 0.8*xi, -2.9*rad, 'Image', /data

;xyouts, -2.*rad, -1.5*rad, 'Lens', /data

return
END 
