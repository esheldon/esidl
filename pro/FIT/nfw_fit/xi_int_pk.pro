pro xi_int_pk,k,pk,r,xi,NumPerPer=NumPerPer,debug=debug,doplot=doplot,$
        xrr=xrr,npd=npd, dospline=dospline
;This works pretty well


nr=n_elements(r)
if nr gt 1 then begin
    xi=dblarr(nr)
    for i=0L, nr-1 do begin
        xi_int_pk,k,pk,r[i],xii,NumPerPer=NumPerPer,npd=npd,dospline=dospline
        xi[i]=xii
    endfor
    return
endif

if n_elements(dospline) eq 0 then dospline=1

if n_elements(NumPerPer) eq 0 then NumPerPer=30.0
if n_elements(npd) eq 0 then npd=100L
;Number of samples per period, should always be > 5 at least

L_wiggle=170d                   ;depends on cosmology but this is ballpark                  
Pref=1d /(2d * !dpi^2 * r)
dk1=2d *!dpi/(NumPerPer*L_wiggle)
dk2=2d *!dpi/(NumPerPer*r)
dktry=min([dk1,dk2])

kmin=0d
kmax=2d
numk=long((kmax-kmin)/dktry)
dk=(kmax-kmin)/double(numk-1)

if keyword_set(debug) then begin
    print,'r= ',r
    if dk1 lt dk2 then print,'Baryon wiggle limited' else print,'Sin(k*r) limited'
    print,'Numk=',numk
    print,'dk= ',dk
endif

kk=dindgen(numk)*dk+kmin

Pkk=interpol(double(Pk),double(k),kk,spline=dospline)
tab=Pkk*kk*sin(kk*r)

symsize=0.5
if keyword_set(doplot) then begin
    pmo=0
    !p.multi=pmo
    w=where(k gt 1.e-8)
    pplot, k[w], Pk[w], psym=8, xrange=[1.e-8,kmax], xstyle=1, /xlog,$
        symsize=symsize
    pplot,/over,kk,Pkk, psym=8, color='red', $
        symsize=symsize
    key=get_kbrd(1)

    !p.multi=[0,1,2]
    pplot,k,Pk*k^2,xr=[kmin,kmax]
    pplot,/over,kk,Pkk*kk^2,color='red',psym=8, $
        symsize=symsize
    pplot,/over,kk,abs(sin(kk*r))*max(Pkk*kk^2)*0.7,color='blue',psym=-8,$
        symsize=symsize
    pplot,kk,tab/max(tab)*max(Pkk*kk^2)*0.7,color='green',psym=-8, $
        symsize=symsize
    !p.multi=0
    key=get_kbrd(1)
endif

int=int_tabulated(kk,tab)
xi_1=Pref*int

if keyword_set(debug) then print,'Integral = ',int,'   Xi(r) = ',xi_1

;now finish with other method
Pref2=Pref/(r^2)

;now need integral \int_xmax^infinity x P(x/r) sin(x)

xmax=r*kmax
num_dec=5
num_per_dec=npd
numx=num_dec*num_per_dec
x=xmax*10.0^(dindgen(numx)/num_per_dec)
Px=interpol(Pk,k,x/r,spline=dospline)
y=x*Px
ly=alog(y)
n=numx

al=(shift(ly,-1)-ly)/(shift(x,-1)-x)
al[n-1]=al[n-2]

Amp=y/exp(al*x)

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
if keyword_set(doplot) then begin
    !p.multi=[0,1,2]
    plot,x,y,xr=xrr,/xst,psym=8,/yno,/xlog,/ylog,$
        symsize=symsize
    pplot,/over,x,Amp*exp(al*x),color='magenta'
    ;oplot,k*r,k*r*Pk,psym=7
    plot,x,total(d,/cum)/integ,/xlog,xtit='x',ytit='Forward',_extra=ex
    ;plot,x,reverse(total(reverse(d),/cum))/integ,/xlog,xtit='x',ytit='Backward',_extra=ex
    !p.multi=0
endif

xi_2=integ*Pref2
xi=xi_1+xi_2

if keyword_set(debug) then begin
    print,'Xi_1=',xi_1
    print,'Xi_2= ',xi_2
    print,'Xi= ',xi
    print,'tab: '
    ;colprint,kk[0:19],tab[0:19],format='(e,e)'
    ;colprint,kk,tab,format='(e,e)'
    ;colprint,lindgen(n_elements(kk)),kk,Pkk,format='(i0,e,e)'
    ;stop
endif

return
end
