PRO der_sigma2xi,r,dersig,r3,xi,C,cov,covxi,nocor=nocor,nocov=nocov,nobinc=nobinc, xi_interior=xi_interior, rrange=rrange, zmean=zmean, omega_m=omega_m

IF n_params() EQ 0 THEN BEGIN
    print,'-syntax  der_sigma2xi,r,dersig,r3,xi,C,cov,covxi,nocor=nocor,nocov=nocov,nobinc=nobinc, xi_interior=xi_interior, zmean=zmean, omega_m=omega_m'
    return
ENDIF

n=n_elements(r)
;rho=8.45075e-2         ;M_sun h^2 pc^-2

; Critical Density
rho=2.78e-1             ;M_sun h^2 pc^-2

; Use omega to get rho at redshift zero
IF n_elements(omega_m) EQ 0 THEN omega_m=.27
rho=rho*omega_m

;; convert to a mean redshift
IF n_elements(zmean) NE 0 THEN rho = rho*(1.0 + zmean)^3 

;; define xi
xi=fltarr(n)

IF (NOT keyword_set(nocov)) AND n_elements(cov) NE 0  THEN BEGIN
    M=fltarr(n,n)
    nocov=0
ENDIF ELSE BEGIN
    nocov=1
ENDELSE

AAlast=0.0
BBlast=0.0
FOR i=0, n-2 DO BEGIN
    FOR j=i, n-2 DO BEGIN
        a=(dersig(j+1)-dersig(j))/(r(j+1)-r(j))
        b=(dersig(j)*r(j+1)-dersig(j+1)*r(j))/(r(j+1)-r(j))
        AA=sqrt(r(j+1)^2 - r(i)^2) - sqrt(r(j)^2 -r(i)^2)
        BB=alog(  (r(j+1)+sqrt(r(j+1)^2-r(i)^2))/(r(j)+sqrt(r(j)^2-r(i)^2)))
        xi(i)=xi(i)-(a*AA+b*BB)/!pi
        IF nocov EQ 0 THEN BEGIN
            M(j,i)=(-AA+BB*r(j+1))/(r(j+1)-r(j))
            IF j NE i THEN M(j,i)=M(j,i)+(AAlast-BBlast*r(j-1))/(r(j)-r(j-1))
            M(j,i)=(-1/!pi)*M(j,i)
            AAlast=AA
            BBlast=BB
        ENDIF
    ENDFOR
    IF nocov EQ 0 THEN BEGIN
                                ;put the last j element in
        M(j,i)=(AAlast-BBlast*r(j-1))/(r(j)-r(j-1))
        M(j,i)=(-1/!pi)*M(j,i)
    ENDIF
ENDFOR

;xi2=M##dersig
;dif=(xi-xi2)

IF nocov EQ 0 THEN BEGIN
    covxi=M##(cov##transpose(M))
ENDIF

IF keyword_set(nocor) THEN BEGIN
    xi=xi/rho
    r3=r(0:n-2)
    xi=xi(0:n-2)
    return
ENDIF

;the rest if for bin correction and C calculation

ymin=1e-6



IF n_elements(rrange) NE 0 THEN BEGIN 
    xmin = rrange[0]
    xmax = rrange[1]
ENDIF ELSE BEGIN 
    w=where(xi GT ymin)
    plot,r(w),xi(w),/xlog,/ylog,psym=1
    print,'click start'
    cursor,xmin,yjunk,/data,/down
    oplot,[1,1]*xmin,[ymin,1000000],color=!blue
    print,'click end'
    cursor,xmax,yjunk,/data,/down
    oplot,[1,1]*xmax,[ymin,1000000],color=!blue
    wait,.2
ENDELSE 
w=where(r GT xmin AND r LT xmax,nfit)
wgood=w
func='power_law_func'

IF nocov THEN begin
    err=replicate(1.0,n)
ENDIF ELSE begin
    err=fltarr(n)
    FOR i=0, n-1 DO err(i)=sqrt(covxi(i,i))
ENDELSE

guess=[1.0,-1.802]
p=mpfitfun(func,r(w),xi(w),err(w),guess,perror=perror,covar=covar,/quiet)
IF n_elements(rrange) EQ 0 THEN oplot,r,power_law_func(r,p),color=!green

rmax=r(n-1)
gam=-p(1)

;use a power series approximation rather than doing
;numerical integration

konst=p(0)*rmax^(-gam)*(gam-1)*gamma((gam-1)/2)/(gamma(1/2.0)*gamma(gam/2))
an=[1.0,1/2.0,3/8.0,5/16.0,35/128.0,63/256.0,231/1024.0,429/2048.0,$
6435/32768.0,12155/65536.0,46189/262144.0,88179/524288.0]

nser=n_elements(an) 
index=2*indgen(nser)
C=fltarr(n)
FOR i=0, n-2 DO BEGIN
    cseries=(r(i)/rmax)^(index) * an /(gam+index)
    
    tsum=total(cseries)
    sumdif=cseries(nser-1)/tsum
    IF sumdif GT .01 THEN begin
        print,'WARNING relative error for partial sums is',sumdif
        print,'not yet converged'
    ENDIF

    C(i)=konst*tsum
ENDFOR

IF n_elements(rrange) EQ 0 THEN BEGIN 
    ;;w=where(xi+C GT ymin)
    ;;plot,r(w),xi(w)+C,/xlog,/ylog,psym=1
    IF 0 THEN begin
        print,'click start'
        cursor,xmin,yjunk,/data,/down
        oplot,[1,1]*xmin,[ymin,1000000],color=!blue
        print,'click end'
        cursor,xmax,yjunk,/data,/down
        oplot,[1,1]*xmax,[ymin,1000000],color=!blue
        w=where(r GT xmin AND r LT xmax,nfit)
    ENDIF
ENDIF 
w=wgood

func='power_law_func'
p=mpfitfun(func,r(w),xi(w)+C(w),err,guess,perror=perror,covar=covar,/quiet)
gam=-p(1)
konst=p(0)*rmax^(-gam)*(gam-1)*gamma((gam-1)/2)/(gamma(1/2.0)*gamma(gam/2))

pl=power_law_func(r,p)
IF n_elements(rrange) EQ 0 THEN BEGIN
    oplot,r,pl,color=!red
    ;;legend,['',''],line=[0,0],color=[!green, !red]
ENDIF 
FOR i=0, n-2 DO BEGIN
    cseries=(r(i)/rmax)^(index) * an /(gam+index)
    C(i)=konst*total(cseries)
ENDFOR

;do last point

C(n-1)=p(0)*r(n-1)^(-gam)
;do bin correction

IF NOT keyword_set(nobinc) THEN BEGIN
    wp=p(0)*(1-gam)*gamma(.5)*gamma((gam-1)/2)/gamma(gam/2) * r^(-gam)
    der_sigma2xi,r,wp,r3out,xiout,cjunk,cov,covjunk,/nobinc, rrange=rrange, zmean=zmean, omega_m=omega_m
    rat=(xiout*rho)/pl(0:n-2)

    IF n_elements(rrange) EQ 0 THEN BEGIN 
        plot,r3out,rat,psym=1,/xlog,tit='Bin correction',ytit='XI_in / XI_out',/yno
        print,'click rmin for bin correction'
        cursor,xmin,yjunk,/data,/down
        oplot,[1,1]*xmin,[ymin,1000000],color=!blue
        print,'click rmax for bin correction'
        cursor,xmax,yjunk,/data,/down
        oplot,[1,1]*xmax,[ymin,1000000],color=!blue
    ENDIF ELSE BEGIN 
        xmin = rrange[0]
        xmax = rrange[1]
    ENDELSE 
    w=where(r GT xmin AND r LT xmax)
    binc=median(rat(w))
    print,'bin correction',binc
    xi=(xi/binc)
    IF nocov EQ 0 THEN covxi=covxi/(binc^2)
ENDIF

C = (C/rho) > 0
xi=xi/rho

;; keep the interior xi as well.  Will use xi_interior-xierr as the
;; lower bound
xi_interior = xi
xi=(xi+C)

IF nocov EQ 0 THEN covxi=covxi/(rho^2)

;now trim off the last element which is undetermined
r3=r(0:n-2)
xi=xi(0:n-2)
xi_interior = xi_interior[0:n-2]
c=c(0:n-2)
IF nocov EQ 0 THEN covxi=covxi(0:n-2,0:n-2)

return
END









