PRO  new_red_fit, cZ, sigma, mag, coeffs, predicted, degree

;*****************************************************************************
;  Procedure:  red_fit
;  
;  Purpose:  Find a first or second degree polynomial fit to redshift data given 
;            as a function of the color magntudes r=X[1],g=X[2],u=X[3],
;            i=X[4],z=X[5]. 
;       
;  Inputs:  cZ:        an array of redshift values of length Num2
;	    sigma:    the standard dev. for each value of cZ
;           mag:        the color magnitudes for each cZ. mag=dblarr[6,Num2]
;           degree:   degree of fit, 1 or 2
;
;  Outputs: coeffs:   an array containing the coefficients for the polynomial
;                     Ask Erin if you need the order
;           predicted:  an array of size Num2 which contains a predicted 
;		        cZ for each of the given redshifts.
;
;  Notes:   The purpose of the array 'indices' may not be apparent.  Please
;           talk to Erin for an explanation
;
;
;  Author:  Erin Sheldon
;
;*****************************************************************************



;*****************************************************************************
;check to see if all the parameters are given
;*****************************************************************************

	IF N_PARAMS() eq 0 THEN BEGIN

		PRINT,"Syntax - red_fit, cZ, sigma, mag, coeffs, predicted, degree"
		return
	END

;*****************************************************************************
;build array Num2, Num2[1] = number of redshifts
;*****************************************************************************

	Num1 = 6              ;number of X's, counting Xo=1.0
	Num2 = size(cZ)

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
		num_coeff = 6
		ENDIF ELSE $
	IF (degree eq 2) THEN BEGIN  
		num_coeff = 21
		ENDIF ELSE  BEGIN

		print,"degree = 1 or 2"
		return
	     ENDELSE


;*****************************************************************************
;prepare a set of indices
;*****************************************************************************
	
        indices = intarr(21,2)
	k=0
	FOR j=0,Num1-1 DO BEGIN
		FOR i=j,5 DO BEGIN
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
	        
		alpha[k]=alpha[k] + cZ[n]*X[n,i]*X[n,j]/(sigma[n]^2)

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

coeffs = TRANSPOSE(alpha) # INVERT(beta)  ;note wierd matrix multiplication


;*****************************************************************************
;make an array of predicted redshifts
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


