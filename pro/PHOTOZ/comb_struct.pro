pro   comb_struct,oldm, oldcz,oldczerr, addm, addcz,addczerr,new_struct

        IF N_PARAMS() eq 0 THEN BEGIN
                
                PRINT,"Syntax - comb_struct, oldm, oldcz, oldczerr, addm, addcz, addczerr, new_struct"
                return
        END


	num1 = n_elements(oldcz)
	num2 = n_elements(addcz)
	tot  = num1+num2

new_struct = create_struct("M", dblarr(6,tot), "CZ", dblarr(tot),"CZERR",dblarr(tot) )

FOR i=0,num1-1 DO BEGIN

	new_struct.m[0,i] = 1.0
	new_struct.m[1,i] = oldm[1,i]
	new_struct.m[2,i] = oldm[2,i]
	new_struct.m[3,i] = oldm[3,i]
	new_struct.m[4,i] = oldm[4,i]
	new_struct.m[5,i] = oldm[5,i]
	
	
ENDFOR

FOR i=0,num1-1 DO BEGIN

	new_struct.cz[i] = oldcz[i]
	new_struct.czerr[i]=oldczerr[i]

ENDFOR

FOR i=0,num2-1 DO BEGIN  

        new_struct.m[0,num1+i] = 1.0   
        new_struct.m[1,num1+i] = addm[1,i]
        new_struct.m[2,num1+i] = addm[2,i]
        new_struct.m[3,num1+i] = addm[3,i]
        new_struct.m[4,num1+i] = addm[4,i]
        new_struct.m[5,num1+i] = addm[5,i]

ENDFOR
	
FOR i=0,num2 -1 DO BEGIN 

        new_struct.cz[num1+i] = addcz[i]
	new_struct.czerr[num1+i] = addczerr[i] 

ENDFOR

	
return
END

	
