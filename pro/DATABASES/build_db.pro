PRO build_db, strin

  str = strin

  str.petrocounts = str.petrocounts - str.reddening
  str.counts_model = str.counts_model - str.reddening

  !priv=2
  dbcreate, 'stripe10_spectra', 1, 1
  dbopen,'stripe10_spectra',1
  
  eq2survey, str.ra, str.dec, lambda, eta

  s=sort(lambda)

  n=n_elements(str)

  cat_no = lindgen(n)

  dbbuild, $
    str[s].run, $
    str[s].rerun, $
    str[s].camcol, $
    str[s].field, $
    str[s].id, $
    str[s].parent, $
    str[s].nchild, $
    str[s].objc_type, $
    str[s].type, $
    str[s].flags, $
    str[s].flags2, $
    str[s].objc_flags, $
    str[s].objc_flags2, $
    str[s].status, $
    str[s].ra/15., $
    str[s].dec, $
    lambda[s], $
    eta[s], $
    str[s].primtarget, $
    str[s].sectarget, $
    str[s].objc_rowc, $
    str[s].objc_colc, $
    str[s].rowc, $
    str[s].colc, $
    str[s].petrocounts, $
    str[s].counts_model, $
    str[s].petror50, $
    str[s].petror90, $
    str[s].z1d, $
    str[s].z1d_error

  return
END 
