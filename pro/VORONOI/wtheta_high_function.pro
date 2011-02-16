FUNCTION wtheta_high_function, x

; g-band
;  a = [0.614798,      2.08658,    -0.570492,  -0.00863386]
; r-band
;  a = [0.654332,      2.20004,    -0.591841,  -0.00831866]
; i-band
;  a = [0.631636,      2.11659,    -0.576719,  -0.00891381]

  ;; this is "gompertz" from comfit
  ;;g = (1 + (x/a(1))^2)
  ;;return, a(0)*(g)^(a(2)/2.0) + a(3)
    
  ;; New stuff. Using power law
  g = (1d + x^2)
  return, !ahigh[0]*(g)^(!ahigh[1]/2d)

END 

