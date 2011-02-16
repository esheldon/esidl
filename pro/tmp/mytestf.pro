FUNCTION mytestf, x, y

  IF n_elements(y) EQ 0 THEN y = 0.
  return, 1./(1.+x^2 + y^2)

END 
