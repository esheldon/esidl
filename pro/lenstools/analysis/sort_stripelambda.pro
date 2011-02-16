PRO sort_stripelambda, lcat, lambda, eta, issouth

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: sort_stripelambda, lcat, lambda, eta[, issouth]'
      return
  ENDIF 

  IF n_elements(issouth) EQ 0 THEN issouth=0

  IF issouth THEN BEGIN
      print
      print,'This is a Southern Stripe.'
      print

      ;; rotate, sort, then rotate back
      rotate_lambda, lambda
      s = sort(lambda)
      lcat = temporary(lcat[s])
      lambda = temporary( lambda[s] )
      eta = temporary( eta[s] )
      ;; rotate back
      rotate_lambda, lambda
  ENDIF ELSE BEGIN 
      ;; sort them
      issouth=0
      s = sort(lambda)
      lcat = temporary(lcat[s])
      lambda = temporary( lambda[s] )
      eta = temporary( eta[s] )
  ENDELSE 

END 
