PRO run_sigmacrit_errors

  ;; you can run sigmacrit_errors_vszl on the outputs

  zLarr = [0.050,0.075,0.100,0.125,0.150,0.175,0.200]

  ;; measure sigmacrit_errors for each zL
  nzL = n_elements(zLarr)

  FOR i=0L, nzL-1 DO BEGIN
      zL = zLarr[i]
      sigmacrit_errors, zL, /dops
  ENDFOR 

END 
