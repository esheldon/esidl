PRO wtheta_get_neighbors_ra, ra, minra, maxra, w, num, issouth=issouth

  ntot=n_elements(ra)
  w=lindgen(ntot)
  
  binary_search, ra, minra, i1, /round
  binary_search, ra, maxra, i2, /round
  
  CASE 1 OF
      (i1 NE -1) AND (i2 NE -1): BEGIN
          w=w[i1:i2]
          num=n_elements(w)
      END 
      (i1 EQ -1) AND (i2 NE -1): BEGIN
          w=w[0:i2]
          num=n_elements(w)
      END 
      (i2 EQ -1) AND (i1 NE -1): BEGIN
          w=w[i1:ntot-1]
          num=n_elements(w)
      END 
      ELSE: BEGIN
          w=-1
          num=0
      END 
  ENDCASE 
  
  return
END 
