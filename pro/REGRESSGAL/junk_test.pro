PRO junk_test,num,ata,btb,av,bv,npair,rsum,clust=clust
COMMON seed,seed
nbins=20L
ata=fltarr(nbins,nbins)    
btb=ata
av=fltarr(nbins) 
bv=fltarr(nbins)
uncert=0.0
npair=fltarr(nbins)
rsum=npair

sigcrit=replicate(1.0,30)
rmax=200.0
rmin=0.0
binsize=(rmax-rmin)/float(nbins)

FOR i=0L, num DO BEGIN
xcen=randomu(seed,3)*500-250.0
ycen=randomu(seed,3)*500-250.0

IF NOT keyword_set(clust) THEN begin
x=randomu(seed,30)*500.0-250.0
y=randomu(seed,30)*500.0-250.0
ENDIF ELSE begin


x=[xcen(0)+50.0*randomn(seed,10),xcen(1)+50.0*randomn(seed,10),xcen(2)+50.0*randomn(seed,10)]
y=[ycen(0)+50.0*randomn(seed,10),ycen(1)+50.0*randomn(seed,10),ycen(2)+50.0*randomn(seed,10)]

ENDELSE

blah1=fltarr(n_elements(x),nbins)
blah2=blah1

IF 0 THEN begin
plot,x,y,psym=1,xr=[-250,250],yr=[-250,250]
oplot,[-rmax,-rmax],[-1000,1000]
oplot,[rmax,rmax],[-1000,1000]
oplot,[-1000,1000],[rmax,rmax]
oplot,[-1000,1000],[-rmax,-rmax]
return
endif
;rad=randomu(seed,30)*200.0
rad=sqrt(x^2+y^2)
;theta=randomu(seed,30)*2*!pi
theta=atan(y,x)
w=where(theta LT 0.0,wif)
IF wif GT 0 THEN theta(w)=2.0*!pi-theta(w)
e1=.8*randomu(seed)-.4
e2=.8*randomu(seed)-.4
;build_system,ata,btb,av,bv,e1,e2,uncert,rad,theta,sigcrit,binsize,rmax,rmin,rsum,npair,nbins,blah1,blah2

build_system,rad,theta,sigcrit,e1,e2,$
      uncert,rmin,rmax,binsize,nbins,ata,btb,av,bv,$
      npair,rsum,blah1,blah2


ENDFOR

return
end
