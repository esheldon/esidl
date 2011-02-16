PRO continued_fraction, in_num, in_denom, coeffs

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: continued_fraction, in_num, in_denom, coeffs'
      return
  ENDIF 

  delvarx, coeffs

  num = in_num
  denom = in_denom
  
  continue=1
  WHILE continue DO BEGIN 
      tcoeff = num/denom
      add_arrval, tcoeff, coeffs
      left = num MOD denom
      IF left EQ 0 THEN return

      num=denom
      denom = left

  ENDWHILE 
      
END 
