PRO  col_mag_data,oldm,oldcol,col,X

	   IF N_PARAMS() eq 0 THEN BEGIN

        	 PRINT,"Syntax - col_mag_data, oldm,oldccol,col,X"
           	 return
           END

	tmpn = n_elements( oldcz )


	only = where( oldcol[*] ne -1 and oldm[1,*] ne -1 and oldm[2,*] ne -1 and $
		      oldm[3,*] ne -1 and oldm[4,*] ne -1 and oldm[5,*] ne -1 )

	num = n_elements(oldcol[only])
	
	X=dblarr(6,num)
	cz=dblarr(num)


	Num = 4  ;Num =4 for 3 variables
	X[0,*] = 1.0
	FOR i=1,Num-1 DO BEGIN
		X[i,*] = oldm[i,only]
	ENDFOR

	col[*] = oldcol[only]

return

END   
	
	


	
