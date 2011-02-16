PRO gal_gal_sim,str,stick=stick
;simulate gal-gal lensing with clustered foregrounds

IF n_params() EQ 0 THEN BEGIN
    print,'-syntax gal_gal_sim,go'
    return
ENDIF
!p.multi=[0,2,2]

IF 1 THEN BEGIN
    ;define
    siz=10000L
    numu=500L
    numc=2L
    numm=10L
    numb=2000L
    wid=100.0
    rref=findgen(siz)
    sig=170.0
    zfore=.1
    zback=.35
    ;dd=200.0
    ;deff=300.0
    dd=angdist_lambda(zfore)
    ds=angdist_lambda(zback)
    dds=angdist_lambda(zback,zfore)
    deff=dd*dds/ds

    tail=10.0
    core=.001
    scale=1.0
    sigcrit=1.0
    uncert=.01
    rmin=50.0
    rmax=2000.0
    nbins=100.0
    binsize=(rmax-rmin)/nbins
    ata=fltarr(nbins,nbins)
    btb=ata
    av=fltarr(nbins)
    bv=av
    npair=lonarr(nbins)
    rsum=av
    buff=100.0
    allow=5000L
    isothermal_distortion,rref,dref,sig,Dd,Deff,tail,Core,scale
    plot,rref,dref
endif
    


xu=randomu(seed,numu)*siz
yu=randomu(seed,numu)*siz

xc=randomu(seed,numc)*siz
yc=randomu(seed,numc)*siz

width=wid*(3.0*randomu(seed,numc)+.2)
mwid=median(width)

xm=-1
ym=-1
FOR i=0, numc-1 DO BEGIN
    nn=numm*((width(i)/mwid)^2)
    nn=round(nn > 5)
    ;print,nn
    xm=[xm,xc(i)+randomn(seed,nn)*width(i)]
    ym=[ym,yc(i)+randomn(seed,nn)*width(i)]

ENDFOR
xm=xm(1:*)
ym=ym(1:*)

xf=[xu,xm]
yf=[yu,ym]
numf=n_elements(xf)
plot,xf,yf,psym=1

xb=randomu(seed,numb)*siz
yb=randomu(seed,numb)*siz

win=where(xb GT (0.0-buff) AND yb GT (0.0-buff) AND $
          xb LT (siz+buff) AND yb LT (siz+buff) )
xb=xb(win)
yb=yb(win)
nb=n_elements(xb)
e1=fltarr(nb)
e2=e1

oplot,xb,yb,psym=7,color=130

close_match,xf,yf,xb,yb,mf,mb,rmax,allow,miss
hh=histogram(mf)
maxmat=max(hh)
print,'mean pair',mean(hh)
nmat=n_elements(mf)
IF maxmat GT allow/2.0 THEN BEGIN
    print,'WARNING maxmatch= ',maxmat,' allow= ',allow
    print,'make ALLOW bigger'
    print,'returning'
    return
endif

um1=uniq(mf,sort(mf))
nun1=n_elements(um1)
um2=uniq(mb,sort(mb))
nun2=n_elements(um2)

FOR i=0L, nun1-1 DO BEGIN
    IF i MOD 100 EQ 0 THEN print,'lensing  ',i,'/',nun1
    w=where(mf EQ mf(um1(i)))
    xhere=xb(mb(w))-xf(mf(um1(i)))
    yhere=yb(mb(w))-yf(mf(um1(i)))

    r=sqrt(xhere^2+yhere^2)
    theta=atan(yhere,xhere)
    wneg=where(theta LT 0.0,wif)
    IF wif GT 0 THEN theta(wneg)=2*!pi+theta(wneg)
    
    isothermal_distortion,r,d,sig,Dd,Deff,tail,Core,scale
    e1(mb(w))=e1(mb(w))-1.0*cos(2.0*theta)*d
    e2(mb(w))=e2(mb(w))-1.0*sin(2.0*theta)*d

    ;tvcircle,2000,xf(mf(um1(i))),yf(mf(um1(i))),/data,color=40

endfor

str={xb:xb/3600.0,yb:yb/3600.0,xf:xf/3600.0,yf:yf/3600.0,e1:e1,e2:e2,$
zb:replicate(zback,numb),zf:replicate(zfore,numf)}


print,'max e1',max(e1)
print,'mean e1',mean(e1)
print,'max e2',max(e2)
print,'mean e2',mean(e2)

plothist,e1,bin=.001
plothist,e2,bin=.001
blah1=fltarr(numf,nbins)
blah2=blah1

FOR i=0L, nun2-1 DO BEGIN
    IF i MOD 100 EQ 0 THEN print,'measuring  ',i,'/',nun2
    w=where(mb EQ mb(um2(i)))
    xhere=xb(mb(um2(i)))-xf(mf(w))
    yhere=yb(mb(um2(i)))-yf(mf(w))
    r=sqrt(xhere^2+yhere^2)
    theta=atan(yhere,xhere)
    wneg=where(theta LT 0.0,wif)
    IF wif GT 0 THEN theta(wneg)=2*!pi+theta(wneg)
    e1err=randomn(seed)*uncert
    e2err=randomn(seed)*uncert

    build_system,r,theta,sigcrit,e1(mb(um2(i)))+e1err,e2(mb(um2(i)))+e2err,$
      uncert,rmin,rmax,binsize,nbins,ata,btb,av,bv,$
      npair,rsum,blah1,blah2

endfor

ainv=invert(ata,status,/double)
IF status NE 0 THEN print,"CAN'T INVERT MATRIX"
binv=invert(btb,status,/double)
IF status NE 0 THEN print,"CAN'T INVERT MATRIX"
a=ainv##transpose(av)
b=binv##transpose(bv)

make_bins,[rmin,rmax],nbins,sta,fin
rr=(sta+fin)/2.0
sh=a+b
cov=ainv+binv
diag=intarr(nbins)
FOR i=0, nbins-1 DO diag(i)=i*(nbins+1)
vars=cov(diag)
err=sqrt(vars)

plot,rr,sh
oplot,rref,dref,color=130

IF keyword_set(stick) THEN ewfewfwe

return
END

