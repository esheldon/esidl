PRO blanton_lumdensity, clr, ld, lderr, omegamat=omegamat

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: blanton_lumdensity, clr, ld, lderr, omegamat=omegamat'
      print,'clr can be any index or array of indices in [0,4]'
      return
  ENDIF 

  IF n_elements(omegamat) EQ 0 THEN omegamat=0.3

  IF omegamat EQ 0.3 THEN BEGIN 
      ld = [4.35, 2.81, 2.58, 3.19, 3.99]*(10d)^(8)
      lderr = [1.08, 0.38, 0.28, 0.37, 0.48]*(10d)^(8)

      ld = ld[clr]
      lderr = lderr[clr]

  ENDIF ELSE BEGIN 
      ld = [4.34, 2.87, 2.74, 3.32, 4.19]*(10d)^(8)
      lderr = [1.10, 0.40, 0.31, 0.40, 0.51]*(10d)^(8)

      ld = ld[clr]
      lderr = lderr[clr]


  ENDELSE 

END 
