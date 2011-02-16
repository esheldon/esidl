FUNCTION wtheta_tot_function, x

  ;; Phil's numbers
  ;a = [0.621, -0.72]
  
  ;; New wtheta stuff in #/Mpc^2
;  a = dblarr(2)
;  a[0] = 0.865306
;  a[1] = -0.671429

  ;; convert to #/kpc^2
;  a[0] = a[0]/(1000d)^2

  g = (1d + x^2)
  return, !atot[0]*(g)^(!atot[1]/2d)
    
END 

