FUNCTION calltot_old, x, y

  r1 = sqrt(x^2 + y^2)
  r2 = sqrt( (x-!conv_x0)^2 + (y-!conv_y0)^2 )
  return, wtheta_tot_FUNCTION(r1)*kappasis_trunc(!conv_sigma,!conv_cutoff,!zsource,!zlens,r2, /core )

END 
