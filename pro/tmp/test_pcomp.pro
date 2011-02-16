PRO test_pcomp, struct, nr, nlum

  nst = n_elements(struct)

  w=where(struct.dn GT 0. AND struct.dn LT 10 AND struct.lum[nlum] GT 0.,nw)


  array = fltarr(3, nw)
  array[0, *] = struct[w].sismass[nr]/1.e14
  array[1, *] = struct[w].lum[nlum]/1.e10
  array[2, *] = struct[w].dn

  coefficients = 1 & eigenvalues = 1 & variances = 1
  result = PCOMP(array, COEFFICIENTS = coefficients, $
                 EIGENVALUES = eigenvalues, VARIANCES =variances)

  PRINT, coefficients
  print
  PRINT, eigenvalues
  print
  PRINT, variances

  print
  print,array[*,0]
  print
  print,result[*,0]

; coeff array converts original data into PCA space
  print
  print,array[*,0] # coefficients

; inverse of coeff array converts variable array back to original data
  print
  print,result[*,0] # INVERT(coefficients)

  return
END 
