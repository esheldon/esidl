PRO ell_parameter_plot

nx=9
ny=9

e1=(findgen(9)-4.)/6.
e2=e1

aplot, 1.0, [-1,1], [-1,1], yrange=[-1,1],xrange=[-1,1],psym=3,$
  ytitle='e!D2',xtitle='e!D1',title='Ellipticity Parameter', $
  xticklen=0.04, yticklen=0.04

oplot,[0,0],[-1,1]
oplot,[-1,1],[0,0]

factor=.1

FOR iy=0L, ny-1 DO BEGIN
    FOR ix=0L, nx-1 DO BEGIN
        
        e=sqrt(e1[ix]^2 + e2[iy]^2)
        IF e LE 0.8 THEN BEGIN
            
            findabtheta, e1[ix], e2[iy], aratio, theta
            
            area = .00225
            a = sqrt(area/aratio)
            
            ;a=.06
            b=aratio*a

            tvellipse, a, b, e1[ix], e2[iy], theta*180./!pi, /data
        ENDIF 
    ENDFOR 
ENDFOR 

return
END 
