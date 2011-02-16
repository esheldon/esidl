FUNCTION wtheta_low_gaussian, x

;; old stuff
; r-band
  ;;a = [-0.0112833,      84.0689,      145.237,  0.000929031]

  ;; new stuff now defined in wtheta_conv
  ;;a = [-5.53803, -0.439464, 0.397702, -0.278953]


  ;; used gaussfit with nterms = 4
  z = (x - !alow[1])/!alow[2]
  return, !alow[0]*exp((-1.0*z^2)/2.0) + !alow[3]

END 
