FUNCTION gmean, r1, r2, dim, double=double

IF n_params() LT 3 THEN BEGIN
    print,'-Syntax: gmean, r1, r2, dim, double=double'
    return,0.
ENDIF 

  e1=dim+1
  e2=dim
  
  IF NOT keyword_set(double) THEN $
    return, float(e2)/e1*( float(r1)^e1 - float(r2)^e1)/ $
                                   (float(r1)^e2 - float(r2)^e2) $
  ELSE $
    return, double(e2)/e1*( double(r1)^e1 - double(r2)^e1)/ $
                                    (double(r1)^e2 - double(r2)^e2)

END 
      
