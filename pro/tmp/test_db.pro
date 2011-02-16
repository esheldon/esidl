PRO test_db, str

  dbopen,'example',1
  
  eq2survey, str.ra, str.dec, lambda, eta

  s=sort(lambda)

  n=n_elements(str)

  cat_no = lindgen(n)

  dbbuild, cat_no, str[s].run, str[s].rerun, str[s].camcol, str[s].field, str[s].id, $
    str[s].ra, str[s].dec, lambda[s], eta[s]

  return
END 
