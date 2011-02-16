;; Flat lambda universe
FUNCTION rhocrit, z=z, h=h, omega_m = omega_m

  IF n_elements(z) EQ 0 THEN z=0.0
  if n_elements(omega_m) eq 0 then omega_m = 0.3
  if n_elements(h) eq 0 then h=1.0

  H2=omega_m *(1.0 + z)^3 + (1-omega_m)
  rhocrit=!rhocrit*H2*h^2

  return,rhocrit

END   
