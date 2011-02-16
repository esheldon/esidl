PRO delta_sigma2xi,r,ds,cov,r3,xi,covxi,xierr, rrange=rrange

IF n_params() EQ 0 THEN BEGIN
    print,'-syntax delta_sigma2xi,r,ds,cov,r3,xi,covxi,xierr, rrange=rrange'
    return
endif

!p.multi=0
print,'Calculating derivative'
delta_sigma2der_sigma,r,ds,dersig,cov,covder,derr
print,'Inverting'
der_sigma2xi,r,dersig,r3,xi,C,covder,covxi, rrange=rrange

!p.multi=[0,2,3]
n=n_elements(r)

dserr=fltarr(n)
dersigerr=fltarr(n)
xierr=fltarr(n-1)
FOR i=0, n-1 DO dserr(i)=sqrt(cov(i,i))
FOR i=0, n-1 DO dersigerr(i)=sqrt(covder(i,i))
FOR i=0, n-2 DO xierr(i)=sqrt(covxi(i,i))
dscorr=fltarr(n,n)
dersigcorr=fltarr(n,n)
xicorr=fltarr(n-1,n-1)
FOR i=0, n-1 DO FOR j=0, n-1 DO dscorr(i,j)=cov(i,j)/sqrt(cov(i,i)*cov(j,j))
FOR i=0, n-1 DO FOR j=0, n-1 DO dersigcorr(i,j)=covder(i,j)/sqrt(covder(i,i)*covder(j,j))
FOR i=0, n-2 DO FOR j=0, n-2 DO xicorr(i,j)=covxi(i,j)/sqrt(covxi(i,i)*covxi(j,j))


xtit='R ( Mpc / h )'
ytit=!csym.delta_cap +' '+ !csym.sigma_cap
yr=[.9*min(ds-dserr),max(ds+dserr)*1.1]
plot,r,ds,psym=1,/xlog,/ylog,xtit=xtit,ytit=ytit,tit=tit,YTICKFORMAT='(F7.1)',yr=yr,/yst
errplot,r,ds-dserr,ds+dserr
tvim2,dscorr,tit='Correlation Matrix'

ytit='- d '+' '+ !csym.sigma_cap +' / dR'
yr=[.9*min(-dersig-dersigerr) > 1e-1,max(-dersig+dersigerr)*1.1]
plot,r,-dersig,psym=1,/xlog,/ylog,xtit=xtit,ytit=ytit,tit=tit,YTICKFORMAT='(F7.1)',yr=yr,/yst
errplot,r,-dersig-dersigerr,-dersig+dersigerr
tvim2,dersigcorr,tit='Correlation Matrix'

xtit='r ( Mpc / h )'
ytit=!csym.xi
yr=[.9*min(xi-xierr),max(xi+xierr)*1.1]
plot,r3,xi,psym=1,/xlog,/ylog,xtit=xtit,ytit=ytit,tit=tit,yr=yr,YTICKFORMAT='(F7.1)',/yst
errplot,r3,xi-xierr,xi+xierr
func='power_law_func'
guess=[1.0,-1.802]
p=mpfitfun(func,r3,xi,xierr,guess,perror=perror,covar=covar,/quiet)
print,p
print,'ro = ',p(0)^(-1.0/p(1)),' +/- ',abs(perror(0)/1.0/p(1))*p(0)^(-1/p(0)-1)
print,'gamma=',-p(1),' +/- ',perror(1)
oplot,r3,power_law_func(r3,p),color=!green
gam=-p(1)
gam=strcompress(gam,/rem)
gam=strmid(gam,0,4)
xyouts,.4,1000,/data,!csym.gamma+" = "+gam,charsize=1.5
tvim2,xicorr,tit='Correlation Matrix'


!p.multi=0

return
end
