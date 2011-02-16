PRO htm_write_neighbors, file, ind1, ind2, hdr=hdr

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: htm_write_neighbors, file, ind1, ind2, hdr=hdr'
      return
  ENDIF 
  
  print
  print,'Outputting file: '+file
  print

  mwrfits, ind1, file, hdr, /create
  mwrfits, ind2, file

  return
END 
