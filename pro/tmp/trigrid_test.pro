PRO tmp, shear, cenx, ceny, kappa, nx=nx, ny=ny


x = shear.dec - cenx
y = shear.ra  - ceny
R = sqrt(x^2 + y^2) > 20./3600. ;degrees
R2 = R^2


kap = -( shear.e1*( x^2 - y^2 ) + shear.e2*2.*x*y )/R^2

triangulate, x, y, tr, b

kappa = trigrid(x, y, kap, tr, nx=nx, ny=ny)

s=120/3600.
makegauss,gauss, [6,6],1.7
use = convol(kappa, gauss, /edge_truncate)
tvim3,use

key=get_kbrd(1)
rotate_plot, use


return

END 
