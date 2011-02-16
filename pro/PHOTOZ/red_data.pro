PRO  red_data,oldm,oldcz,oldczerr,X,cz,czerr

	   IF N_PARAMS() eq 0 THEN BEGIN

        	 PRINT,"Syntax - red_data, oldm, oldcz, oldczerr, X, cz, czerr"
           	return
           END

	tmpn = n_elements( oldcz )


	only = where( oldcz[*] ne -1 and oldm[1,*] ne -1 and oldm[2,*] ne -1 and $
		      oldm[3,*] ne -1 and oldm[4,*] ne -1 and oldm[5,*] ne -1 )

	num = n_elements(oldcz[only])
	
	X=dblarr(6,num)
	cz=dblarr(num)
	czerr = dblarr(num)

	X[0,*] = 1.0
	FOR i=1,5 DO BEGIN
		X[i,*] = oldm[i,only]
	ENDFOR

	cz[*] = oldcz[only]/3.0e5
	czerr[*] = oldczerr[only]/3.0e5

return

END   
	
	


	
