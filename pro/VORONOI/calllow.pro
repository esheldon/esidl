FUNCTION calllow, x, y

  r1 = sqrt(x^2 + y^2)
  r2 = sqrt( (x-!conv_x0)^2 + (y-!conv_y0)^2 )
  return, wtheta_low_gaussian(r1)*sigmasis_trunc(!conv_sigma,!conv_cutoff,r2, /core)

END 
