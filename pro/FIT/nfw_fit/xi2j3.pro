pro xi2j3,rr,xi,j3
;integrate a tabulated xi to get J3
;J3(r)*rho_bar = M(r)
;ignore the last point which is meaningless
;Check this!!!!!!!!

if n_params() eq 0 then begin
    print,'-syntax xi2j3,r,xi,j3'
    return
endif

num=n_elements(xi)
nnr=n_elements(rr)
r=rr
if nnr eq num+1 then r=r[0:num-1]

wneg=where(xi lt 0,nneg)

if nneg gt 0 then begin
    print,'Error - xi must be positive for this powerlaw interpolation'
    return
endif

lr=alog(r)
lx=alog(xi)

al=(lx-shift(lx,1))/(lr-shift(lr,1))
al[0]=al[1]

A=xi/r^(al)

ex=3.0+al
Rb=r
Ra=shift(r,1)
Ra[0]=0.0
int0=A*(1.0/ex) *(Rb^ex -Ra^ex)

j3=4*!pi*total(int0,/cum) 

return
end
