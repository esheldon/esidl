PRO  truered_fit, logz, sigma, degree, delta, coeffs 

;*****************************************************************************
;check to see if all the parameters are given
;*****************************************************************************

	IF N_PARAMS() eq 0 THEN BEGIN

		PRINT,"Syntax - red_fit, logz, sigma, degree,delta,coeffs"
		return
	END

;*****************************************************************************
;build array N
;*****************************************************************************

	s = size(logz)
	N = intarr(6)

	FOR i=1,5 DO  BEGIN   

		N[i] = s[i]
	ENDFOR

	N[0]=0
	N[0]=MAX(N)

;*****************************************************************************
;build array X, the color magnitude grid
;*****************************************************************************

	X = DBLARR(6,N[0])

	FOR i=1,5 DO  BEGIN 

		X[i,*] = 1.0*FINDGEN( N[0] )
	ENDFOR

	FOR i=0,N[0]-1  DO BEGIN 

		X[0,i] = 1.0
	ENDFOR

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
	FOR j=0,5 DO BEGIN
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
	a     = DBLARR(num_coeff)
	M     = INTARR(num_coeff)
	M[0]  = 0
;*****************************************************************************
;calculate the vector alpha
;*****************************************************************************

	FOR k=0,num_coeff-1 DO BEGIN
	  
	  i = indices[k,0] 
	  j = indices[k,1]

	   FOR i1=0, N[1]-1 DO BEGIN
	      M[1]=i1	 
	      FOR i2=0, N[2]-1 DO BEGIN
                 M[2]=i2
                 FOR i3=0, N[3]-1 DO BEGIN
                    M[3]=i3
                    FOR i4=0, N[4]-1 DO BEGIN
                       M[4]=i4
                       FOR i5=0, N[5]-1 DO BEGIN
                          
                          M[5]=i5
                          tmp = logz[ M[1],M[2],M[3],M[4],M[5] ]
			  tmp = tmp*X[i,M[i] ] * delta[i] *X[j,M[j] ] *delta[j]	
                          tmp = tmp/sigma[ M[1],M[2],M[3],M[4],M[5] ]
                          alpha[k] = alpha[k] + tmp

                       ENDFOR
                    ENDFOR
                 ENDFOR
              ENDFOR
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

              FOR i1=0, N[1]-1 DO BEGIN
                 M[1]=i1
                 FOR i2=0, N[2]-1 DO BEGIN
                    M[2]=i2
                    FOR i3=0, N[3]-1 DO BEGIN
                       M[3]=i3
                       FOR i4=0, N[4]-1 DO BEGIN
                          M[4]=i4
                          FOR i5=0, N[5]-1 DO BEGIN

                             M[5]=i5
                             tmp = X[i, M[i] ]*X[j, M[j] ]*delta[i]*delta[j]
                             tmp = tmp*X[u, M[u] ]*X[v, M[v] ]*delta[u]*delta[v]
                             tmp = tmp/sigma[ M[1],M[2],M[3],M[4],M[5] ]
                             beta[k,l] = beta[k,l] + tmp

                          ENDFOR
                       ENDFOR
                    ENDFOR
                 ENDFOR
              ENDFOR
	         beta[l,k] = beta[k,l]		
           ENDFOR
        ENDFOR

coeffs = TRANSPOSE(alpha) # INVERT(beta)

return
END
