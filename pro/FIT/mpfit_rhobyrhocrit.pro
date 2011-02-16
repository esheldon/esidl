;; Fit <rho(<r)>/rhocrit for r200, and a power
function mpfit_rhobyrhocrit, x, p

  ;; p[0] = r200
  ;; p[1] = power

  return, 200d*(x/p[0])^p[1]

end 
