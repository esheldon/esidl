FUNCTION prange, v1, v2, v1err, v2err, noerror=noerror, $
                 slack=slack, symmetric=symmetric

  np = n_params()
  IF (np EQ 0) OR (np EQ 3) THEN BEGIN 
      print,'-Syntax: result = prange(v1, [v2,] v1err, [v2err,] [/noerror,/slack,/symmetric])'
      print,' either send (v1,v1err)  or (v1,v2) or all'
      return,-1
  ENDIF 

  IF np EQ 1 THEN BEGIN
      return,[min(v1), max(v1)]
  ENDIF 

  wmax1 = ( where( v1 EQ max(v1)) )[0] 
  wmin1 = ( where( v1 EQ min(v1)) )[0]
  wmax2 = ( where( v2 EQ max(v2)) )[0]
  wmin2 = ( where( v2 EQ min(v2)) )[0]

  IF np EQ 2 THEN BEGIN  
      IF NOT keyword_set(noerror) THEN BEGIN ; Assume v2=v1err
          max = v1[wmax1] + v2[wmax1]
          min = v1[wmin1] - v2[wmin1]
          range = [ min, max ]
      ENDIF ELSE BEGIN          ;v1 and v2 seperate vectors  no error bars
          max1 = v1[wmax1]
          min1 = v1[wmin1]
          max2 = v2[wmax2]
          min2 = v2[wmin2]

          range = [ min([min1,min2]), max([max1,max2]) ]
      ENDELSE 
  ENDIF ELSE BEGIN              ;All must have been sent
      
      max1 = v1[wmax1] + v1err[wmax1]
      min1 = v1[wmin1] - v1err[wmin1]
      max2 = v2[wmax2] + v2err[wmax2]
      min2 = v2[wmin2] - v2err[wmin2]

      range = [ min([min1,min2]), max([max1,max2]) ]
  ENDELSE 

  IF n_elements(slack) NE 0 THEN BEGIN 
      IF slack NE 1 THEN BEGIN 
          fac1=1.0/slack
          fac2=slack
      ENDIF ELSE BEGIN 
          fac1=0.9
          fac2=1.1
      ENDELSE 
      IF range[0] GT 0.0 THEN range[0] = fac1*range[0] ELSE range[0]=fac2*range[0]
      IF range[1] GT 0.0 THEN range[1] = fac2*range[1] ELSE range[1]=fac1*range[1]
  ENDIF 

  ;; if requested make symmetric across zero
  IF keyword_set(symmetric) THEN BEGIN 
      IF (range[0] LT 0.0) AND (range[1] GT 0.0) THEN BEGIN 
          maxy=max(abs(range))
          range=[-maxy, maxy]
      ENDIF ELSE BEGIN 
          tmax = max( abs(range) )
          range = [-tmax, tmax]
      ENDELSE 
  ENDIF 

  return, range

END 
