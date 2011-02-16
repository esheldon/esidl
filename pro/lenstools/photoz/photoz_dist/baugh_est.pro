PRO baugh_est, zc

  ;; the baugh est. thing

  minz = 0.0
  maxz = 2.0

  nz = 100
  z = arrscl( findgen(nz), minz, maxz )

  zmed = 1.412*zc

  zdist = 1.5*z^2/zc^3*exp(-(z/zc)^1.5)

  plot,z,zdist
  oplot, [zmed, zmed], [0,1000], color=!blue
  oplot, [zc, zc], [0,1000], color=!red

  npts = 100
  gauleg, minz, maxz, npts, ZZi, WWi

  zdist_vals = interpol(zdist, z, ZZi)
  norm = total( zdist_vals*WWi )

  print,norm

END 
