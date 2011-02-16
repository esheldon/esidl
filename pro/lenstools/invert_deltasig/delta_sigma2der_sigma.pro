PRO delta_sigma2der_sigma,r,delsig,dersig,covdel,covder,err,plot=plot,mat=mat
;calculate the derivative of sigma from
;finite difference equation of delta sigma
;also transform covariance matrix 
;results  array rr and dersig that are 2 elements smaller 
;since it computes the symetric derivative
;and throws away the two endpoints
;also trims off teh edges of covder
;inputs:
;r         the radius (any units)
;delsig    the delta sigma measurements
;covdel    the covariance matrix of delsig
;outputs:   
;rr        the outputs radius that has been trimmed on each end
;dersig    the derivitive of sigma          
;covder    the covariance matrix of dersig propagated from covdel
;err       the sqrt(diagonal) of covderr
;keywords:
;plot      will plot both r,delsig and rr,dersig with error bars
  
IF n_params() EQ 0 THEN BEGIN
    print,'-syntax delta_sigma2der_sigma,r,delsig,dersig,covdel,covder,err,plot=plot,mat=mat'
    return
endif

n=n_elements(r)
sdel=shift(delsig,-1)-shift(delsig,1)
sr=shift(r,-1)-shift(r,1)

sdel(0)=delsig(1)-delsig(0)
sr(0)=r(1)-r(0)

sdel(n-1)=delsig(n-1)-delsig(n-2)
sr(n-1)=r(n-1)-r(n-2)

dersig=-2*delsig/r - sdel/sr

mat=fltarr(n,n)
;the transformation matrix

FOR i=0, n-1 DO BEGIN
    mat(i,i)=-2/r(i)
    IF i LT (n-1) then mat(i+1,i)=-1/sr(i)
    IF i GT 0 THEN mat(i-1,i)=1/sr(i)
    IF i EQ 0 THEN BEGIN
        mat(0,0)= -2/(r(0)) + 1/(sr(0))
        mat(1,0)= -1/sr(0)
    ENDIF
    IF i EQ n-1 THEN BEGIN
        mat(n-1,n-1)= -2/r(n-1) -1/sr(n-1)
        mat(n-2,n-1)= 1/sr(n-1)
    ENDIF
ENDFOR

IF n_elements(covdel) GT 0 THEN begin
    covder=mat##(covdel##transpose(mat))
endif

check=mat##delsig
dif=abs((dersig-check)/dersig)
IF max(dif) GT .001 THEN BEGIN
    print,'WARNING matrix formulation gives different answer'
    print,'check the code'
    print,max(dif)
ENDIF


!p.multi=[0,1,2]

IF n_elements(covdel) GT 0 THEN begin
    err=fltarr(n)
    FOR i=0, n-1 DO BEGIN
        err(i)=sqrt(covder(i,i))
    ENDFOR
ENDIF

IF keyword_set(plot) THEN BEGIN
    xtit="R"
    ytit="Delta Sigma"
    delerr=fltarr(n)
    FOR i=0, n-1 DO delerr(i)=sqrt(covdel(i,i))
    !p.multi=[0,1,2]
    plot,r,delsig,psym=1,/xlog,/ylog,xtit=xtit,ytit=ytit
    ytit=" - d Sigma / dR "
    errplot,r,delsig-delerr,delsig+delerr
    plot,r,-dersig,psym=1,/xlog,/ylog,xtit=xtit,ytit=ytit
    errplot,r,-dersig-err,-dersig+err
ENDIF

return
end
