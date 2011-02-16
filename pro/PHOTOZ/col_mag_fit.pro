PRO  col_mag_fit, col, sigma,mag,coeffs,predicted, degree

;*****************************************************************************
;  Procedure:  col_mag_fit
;  
;  Purpose:  Find a first or second degree polynomial fit to redshift data given 
;            as a function of the color magntudes r=X[1],g=X[2],u=X[3],
;            i=X[4],z=X[5]
;       
;  Inputs:  col:      an array of color magnitudes of length Num2
;	    sigma:    the standard dev. for each value of cZ
;           mag:        the color magnitudes for each col. NOTE X[0]=1.0
;           degree:   degree of fit, 1 or 2
;
;  Outputs: coeffs:   an array containing the coefficients for the polynomial
;                     Ask Erin if you need the order
;           predicted:  an array of size Num2 which contains a predicted 
;		        col for each of the given ones.
;
;  Author:  Erin Sheldon
;
;*****************************************************************************



;*****************************************************************************
;check to see if all the parameters are given
;*****************************************************************************

	IF N_PARAMS() eq 0 THEN BEGIN

		PRINT,"Syntax - col_mag_fit, col, sigma, mag, coeffs, predicted, degree"
		return
	END

;*****************************************************************************
;build array N, N[1] = number of redshifts
;*****************************************************************************

	Num1 = 4             ;number of X's, counting Xo
	Num2 = size(col)

;*****************************************************************************
;The array mag needs to be put into the correct format for the rest of the
;code.  I wrote the code for an array of this type
;*****************************************************************************


  IF (Num2[1] ge Num1) THEN BEGIN

        X=dblarr(Num2[1],Num1)

        X[*,0]=1.0
        FOR i=1,Num1-1 DO BEGIN

                X[*,i]=mag[i,*]
        ENDFOR

  END

;*****************************************************************************
;determine how many coefficients are needed in the polynomial
;*****************************************************************************

	IF (degree eq 1) THEN BEGIN
		num_coeff = 4
		ENDIF ELSE $
	IF (degree eq 2) THEN BEGIN  
		num_coeff = 10
		ENDIF ELSE  BEGIN

		print,"degree = 1 or 2"
		return
	     ENDELSE


;*****************************************************************************
;prepare a set of indices
;*****************************************************************************
	
        indices = intarr(10,2) ;we need all 10
	k=0
	FOR j=0, Num1-1  DO BEGIN
		FOR i=j, Num1-1  DO BEGIN
			indices[k,0]=j
			indices[k,1]=i
                        k=k+1
		ENDFOR

	ENDFOR

;*****************************************************************************
;some declarations, see header      dbl arrays initialized to zero
;*****************************************************************************

	alpha = DBLARR(num_coeff)
	beta  = DBLARR(num_coeff,num_coeff)

;*****************************************************************************
;calculate the vector alpha
;*****************************************************************************

	FOR k=0,num_coeff-1 DO BEGIN
	  
	    i = indices[k,0] 
	    j = indices[k,1]
	    FOR n=0, Num2[1]-1 DO BEGIN
	        
		alpha[k]=alpha[k] + col[n]*X[n,i]*X[n,j]/(sigma[n]^2)

	    ENDFOR
        ENDFOR

;*****************************************************************************
;calculate array beta
;*****************************************************************************


	FOR k=0,num_coeff-1  DO BEGIN
           FOR l=k, num_coeff-1 DO BEGIN

	      i = indices[k,0]
              j = indices[k,1]
              u = indices[l,0]
	      v = indices[l,1]

		FOR n=0, Num2[1]-1 DO BEGIN
		   
		    beta[k,l] = beta[k,l] + X[n,i]*X[n,j]*X[n,u]*X[n,v]/(sigma[n]^2)
		
		ENDFOR
	    beta[l,k]= beta[k,l]
	    ENDFOR
        ENDFOR

coeffs = TRANSPOSE(alpha) # INVERT(beta)


;*****************************************************************************
;make an array of predicted color magnitudes
;*****************************************************************************

predicted = dblarr(Num2[1])


	FOR n=0,Num2[1]-1 DO BEGIN
		
		FOR k=0,num_coeff-1 DO BEGIN

			i=indices[k,0]
			j=indices[k,1]
			tmp = coeffs[k]*X[n,i]*X[n,j]
			predicted[n]=predicted[n]+tmp 

		ENDFOR
	ENDFOR


print,coeffs

return
END


