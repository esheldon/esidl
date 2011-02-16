PRO  red_fake,cz,X,num,max_X


 IF N_PARAMS() eq 0 THEN BEGIN

         PRINT,"Syntax - red_fake,cz,X,num,max_X"
           return
           END

a=1.0
b=.03
c=0.13
d=.01
e=.1
f=.05
g=.05
h=.15

nn=6
cz = dblarr(num)
X=dblarr(num,nn)
seed = 2.34

X[*,0]=1.0
FOR i=0,num-1 DO BEGIN
	FOR j=1,nn-1 DO BEGIN
		X[i,j]= randomu(seed)*max_X
	ENDFOR
ENDFOR
      
FOR i=0,num-1 DO BEGIN
					                
            cz[i] = a + b*X[i,1] + c*X[i,2] + d*X[i,3] + e*X[i,4] + f*X[i,5] +$
		    g*X[i,1]^2 + h*X[i,1]*X[i,2] + 2*randomn(seed) 
END 

print,max(cz)

return

END 


	


	
