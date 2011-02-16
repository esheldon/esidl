PRO chisq_diff, v1, v1err, v2, v2err, chisq, redchisq

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; Calculate the chisq difference between two sets of data 
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_params() EQ 0 THEN BEGIN
      print,'-Syntax: chisq_diff, v1, v1err, v2, v2err, chisq, redchisq'
      return
  ENDIF 
      

  n = n_elements(v1)

  chisq = 0.
  usen = 0
  FOR i=0, n-1 DO BEGIN
      
      IF (v1err[i] NE 0) AND (v2err[i] NE 0) THEN BEGIN 
          sig2 = v1err[i]^2 + v2err[i]^2

          chisq = chisq + (v1[i] - v2[i])^2/sig2
          usen = usen+1
      ENDIF 

  ENDFOR 

  redchisq = chisq/usen

  print,'Chisq: ',ntostr(chisq)
  print,'Reduced Chisq: ',ntostr(redchisq)
  print,ntostr(usen),' degrees of freedom (out of ',ntostr(n),')'
return
END 
