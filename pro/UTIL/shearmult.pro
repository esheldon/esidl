PRO shearmult, e1a, e2a, e1b, e2b, e1out, e2out

  ;; shear composition operator

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: shearmult, e1a, e2a, e1b, e2b, e1out, e2out'
      return
  ENDIF 

  na = n_elements(e1a)
  na2 = n_elements(e2a)
  IF na NE na2 THEN message,'All arrays must be the same size'

  nb = n_elements(e1b)
  nb2 = n_elements(e2b)
  IF nb NE nb2 THEN message,'All arrays must be the same size'

  IF na NE nb THEN message,'All arrays must be the same size'

  dotp = e1a*e1b + e2a*e2b
  factor = (1.-sqrt(1-e1b*e1b-e2b*e2b)) / (e1b*e1b + e2b*e2b)

  e1out = (e1a + e1b + e2b*factor*(e2a*e1b - e1a*e2b))/(1+dotp) 
  e2out = (e2a + e2b + e1b*factor*(e1a*e2b - e2a*e1b))/(1+dotp) 

END 
