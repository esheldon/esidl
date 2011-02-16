pro xi_integral,k,Pk,r,xi,d,plot=plot,_extra=ex
;k in h/Mpc
;x is k/r
;pk is power spectrum, must be positive
;r is a scalar in Mpc/h

x=k*r
y=x*Pk
ly=alog(y)
n=n_elements(k)

al=(shift(ly,-1)-ly)/(shift(x,-1)-x)
al[n-1]=al[n-2]

A=y/exp(al*x)

;integral boundaries
a=x
b=shift(x,-1)

norm=y/(1+al^2)
Ta=al*sin(a)-cos(a)
Tb=al*sin(b)-cos(b)
dif=exp(al*(b-a))*Tb-Ta
d=norm*dif
d=d[0:n-2]

integ=total(d)
if keyword_set(plot) then begin
    pmo=!p.multi
    !p.multi=[0,1,2]
    plot,k,total(d,/cum)/integ,/xlog,xtit='k',ytit='Forward',_extra=ex
    plot,k,reverse(total(reverse(d),/cum))/integ,/xlog,xtit='k',ytit='Backward',_extra=ex
    !p.multi=pmo
endif

xi=integ/(2*!pi^2 *r^3)

return
end
