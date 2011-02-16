FUNCTION prob_chisq, degfree, chisq

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: prob = prob_chisq(degfree, chisq)'
      return, -1.0
  ENDIF 

  ;; return the probability that the chisq would exceed the input
  ;; chisq *by chance* even for a correct model

  return, 1. - igamma( degfree/2.0, chisq/2.0 )

END 
