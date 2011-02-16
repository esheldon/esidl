FUNCTION interp_meanscinv, zl, zs, zserr, scinv_struct

  zli = interpol(scinv_struct.zli, scinv_struct.zl, zl)
  zsi = interpol(scinv_struct.zsi, scinv_struct.zs, zs)
  zserri = interpol(scinv_struct.zserri, scinv_struct.zserr, zserr)

  print,zli,zsi,zserri

  return, interpolate(scinv_struct.mean_scinv, $
                      zli, zsi, zserri)

END 
