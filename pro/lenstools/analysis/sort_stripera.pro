PRO sort_stripera, lcat, lra, ldec, issouth

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: sort_stripelambda, lcat, ra, dec, issouth'
      return
  ENDIF 

  IF issouth THEN BEGIN
      print
      print,'This is a Southern Stripe.'
      print

      ;; rotate, sort, then rotate back
      rotate_ra, lra
      s=sort(lra)
     
      lcat = temporary(lcat[s])
      lra = temporary( lra[s] )
      ldec = temporary( ldec[s] )
      ;; rotate back
      rotate_ra, lra
  ENDIF ELSE BEGIN 
      ;; sort them
      issouth=0
      s = sort(lra)
      lcat = temporary( lcat[s] )
      lra = temporary( lra[s] )
      ldec = temporary( ldec[s] )
  ENDELSE 

END 
