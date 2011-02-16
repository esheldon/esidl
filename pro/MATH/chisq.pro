PRO chisq, v1, v1err, v2, chisq

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; Calculate the chisq between data with error (v1) and model (v2)
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_params() EQ 0 THEN BEGIN
      print,'-Syntax: chisq, v1, v1err, v2, chisq'
      return
  ENDIF 
      

  n = n_elements(v1)

  chisq = 0.
  usen = 0
  FOR i=0, n-1 DO BEGIN
      
      IF (v1err[i] NE 0) THEN BEGIN 
          chisq = chisq + (v1[i] - v2[i])^2/v1err[i]^2
          usen = usen+1
      ENDIF 

  ENDFOR 

  print,'Chisq: ',ntostr(chisq)

  print,'Used ',ntostr(usen),' points (out of ',ntostr(n),')'
return
END 
