FUNCTION tmean, array, nelements

  IF nelements GE 1 THEN return, total(array)/nelements ELSE return, array[0]

END 

