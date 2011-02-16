PRO  col_mag_fake,col,X,num,delta


 IF N_PARAMS() eq 0 THEN BEGIN

         PRINT,"Syntax - col_mag_fake,col,X,num,delta"
           return
           END

a=1.0
b=2.0
c=0.5
d=1.0


Num3 = 4  ;4 for 3 variables
col = dblarr(num)
max = num*delta
X=dblarr(num,Num3)
seed = 3.218


X[*,0]=1.0
FOR i=0,num-1 DO BEGIN
	FOR j=1,Num3-1 DO BEGIN
		X[i,j]= randomu(seed)*max
	ENDFOR
ENDFOR
      
FOR i=0,num-1 DO BEGIN
					                
            col[i] = a + b*X[i,1] + c*X[i,2] + d*X[i,3] 
END

return

END 


	


	
