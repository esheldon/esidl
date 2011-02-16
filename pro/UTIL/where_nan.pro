FUNCTION where_nan, array, count, complement=complement, ncomplement=ncomplement, l64=l64

  IF n_params() LT 1 THEN BEGIN 
      print,'result = where_nan(array, count, complement=complement, ncomplement=ncomplement, l64=l64)'
      return,-1
  ENDIF 

  return, where( finite(array,/nan), count, $
                 complement=complement, ncomplement=ncomplement, $
                 l64=l64 )

END 
