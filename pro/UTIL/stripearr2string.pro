FUNCTION stripearr2string, stripes

  nstripe = n_elements(stripes)

  stripe_string = ''
  FOR ist=0L, nstripe-1 DO BEGIN 
      IF ist EQ nstripe-1 THEN BEGIN 
          stripe_string = stripe_string + stripe2string(stripes[ist])
      ENDIF ELSE BEGIN 
          stripe_string = stripe_string + stripe2string(stripes[ist]) + '_'
      ENDELSE 
  ENDFOR 

  return,stripe_string

END 
