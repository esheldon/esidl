PRO tmp, uncert

n=50000L
dir='/sdss3/usrdevel/esheldon/idl.lib/SHEARSIM/'
h = mrdfits(dir+'e_uncert_hist.fit',1)

IF n_elements(uncert) EQ 0 THEN genrand, h.hist, h.x, n, uncert,/double

fac = ( lindgen(n) MOD 2 )*2 - 1
a = randomu(seed,n)
s=sort(a)
fac = fac[s]

ee=abs( .32*randomn(seed,n)+.55 )
theta = arrscl(a, -!pi/2., !pi/2., arrmin=0., arrmax=1.)
e1 = ee*cos(2*theta)
e2 = ee*sin(2*theta)

plothist, e1, bin=.01

key=get_kbrd(1)

plothist, e1+fac*uncert, bin=.01, /overplot,linestyle=3

return

END 
