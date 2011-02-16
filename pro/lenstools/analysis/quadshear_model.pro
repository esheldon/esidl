pro quadshear_model,ellip,rmax,bin,grid,r,etan,q,$
                 etan1,etan2,erad1,erad2,$
                 theta=theta

;;  Creates a square grid of points, calculates e1,e2 at every point.  
;;Calculate tan and rad shear at every point.  Bins the points in radial
;;bins, Then separates in 2 half spaces, |y|>|x| and |y|<|x|.  Then avg
;;the shear in the bins and the two spaces.
;;
;; Also calculates quadrupole moment E.S.S. 

; theta=theta  The angle by which the galaxy is misaligned.


IF N_PARAMS() eq 0 THEN BEGIN

		PRINT,"Syntax- tanshear_quad,ellip,rmax,bin,grid,r,"
		PRINT,"etan1,etan2,erad1,erad2,theta=theta"        
		return
	END
;;;;
; bin - The size of the radial rings
; grid - The spacing between the grid points
;
;;;;

;; Convert from ellip to e
rat=sqrt( (1-ellip)/(1+ellip) )
e= (1-rat)/(1+rat)
	

num=rmax/grid*2.0
numrad=rmax/bin
if (rmax mod bin ne 0) then numrad=numrad+1
x=(findgen(num)+.5)*grid-rmax
y=(findgen(num)+.5)*grid-rmax
r=(indgen(numrad)+1)*bin*1.0

etan = fltarr(numrad)
q=etan
qr=q
qi=q
qradr=q
qradi=q

etan1=fltarr(numrad)
etan2=fltarr(numrad)
erad1=fltarr(numrad)
erad2=fltarr(numrad)
sum1=intarr(numrad)
sum2=intarr(numrad)

;;debug
etanall=fltarr(num,num)

IF Keyword_set(theta) THEN a=theta ELSE a=0

FOR i=0,num-1 DO BEGIN
   FOR j=0, num-1 DO BEGIN
     IF (x(i)^2+y(j)^2 ge 10) THEN BEGIN

; Account for misalignement
	
	xi=x(i)*cos(-a)+y(j)*sin(-a)
  	yj=y(j)*cos(-a)-x(i)*sin(-a)
	
	siemd_shear,e,xi,yj,e1,e2

	ne1=e1*cos(2*a)+e2*sin(2*a)
	ne2=e2*cos(2*a)-e1*sin(2*a)
	e1=ne1
	e2=ne2

	etanall(i,j)=e1*(y(j)^2-x(i)^2)/(x(i)^2+y(j)^2) - 2*e2*x(i)*y(j)/(x(i)^2+y(j)^2)
	erad=e1*2*x(i)*y(j)/(x(i)^2+y(j)^2) +e2*(y(j)^2-x(i)^2)/(x(i)^2+y(j)^2)
	radius=sqrt(x(i)^2+y(j)^2)
;; Find in which ring this point falls, then on which quadrant
	IF (radius le rmax) THEN BEGIN
	  k=0
	  WHILE radius gt r(k) do k=k+1 
	  IF ( abs(y(j)) gt abs(x(i)) ) THEN BEGIN
	    etan1(k)=etan1(k)+etanall(i,j)
	    erad1(k)=erad1(k)+erad 
            sum1(k)=sum1(k)+1	;; Keep track of how many in each ring, for mean
	  ENDIF 
	  IF ( abs(y(j)) lt abs(x(i)) ) THEN BEGIN
	    etan2(k)=etan2(k)+etanall(i,j)
	    erad2(k)=erad2(k)+erad 
            sum2(k)=sum2(k)+1
	  ENDIF 
          etan[k] = etan[k] + etanall[i,j]

          ;; now quadrupole terms
          r2= radius^2
          qr[k] = qr[k] + (x[i]^2-y[j]^2)/r2*etanall[i,j]
          qi[k] = qi[k] + 2.0*x[i]*y[j]/r2*etanall[i,j]
          qradr[k] = qradr[k] - 2.0*x[i]*y[j]/r2*erad
          qradi[k] = qradi[k] + (x[i]^2-y[j]^2)/r2*erad

	ENDIF

     ENDIF
   ENDFOR
ENDFOR

sum = sum1+sum2
etan = etan/sum
qr = qr/sum
qi = qi/sum
qradr = qradr/sum
qradi = qradi/sum
;q = sqrt( (qr+qradr)^2 + (qi+qradi)^2 )
q = qr+qradr

;forprint,qr,qradr,qi,qradi

etan1=etan1/sum1
etan2=etan2/sum2
erad1=erad1/sum1
erad2=erad2/sum2

Return
End
