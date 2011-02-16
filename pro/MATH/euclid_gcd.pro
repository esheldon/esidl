FUNCTION euclid_gcd, a, b

  IF b EQ 0 THEN return, a $
    ELSE BEGIN 
      IF a LT b THEN message,'first input must be greater than second'
      return,euclid_gcd(b, a MOD b)
  ENDELSE 

END 
