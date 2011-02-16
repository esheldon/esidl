FUNCTION diagonal_array, diag_elements

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: retval = diagonal_array(diag_elements)'
      return,-1
  ENDIF 

  nn = n_elements(diag_elements)
  tmp = [diag_elements[0]]
  tmp[0] = 0

  retval = replicate(tmp[0], nn, nn)
  
  FOR i=0L, nn-1 DO retval[i,i] = diag_elements[i]

  return, retval

END 
