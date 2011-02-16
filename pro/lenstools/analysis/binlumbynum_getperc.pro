PRO binlumbynum_getperc, nbin, perc

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: binlumbynum_getperc, nbin, perc'
      return
  ENDIF 

  CASE nbin OF
      2: perc =  [0.8, 0.2]
      3: perc = [0.868082, 0.0882179, 0.0436542]
      4: perc = [0.693466, 0.174617, 0.0882179, 0.0436542]
      ELSE: message,"Don't know about nbin = "+ntostr(nbin)
  ENDCASE 

END 
