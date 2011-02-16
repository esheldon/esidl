FUNCTION sortby2, param1, param2

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: s = sortby2(param1, param2)'
      print
      print,' param1,2 must be long/int and same size'
      return
  ENDIF 

  IF (datatype(param1) NE 'INT') AND (datatype(param1) NE 'LON') THEN message,'param 1 is not of type int or long'
  IF (datatype(param2) NE 'INT') AND (datatype(param2) NE 'LON') THEN message,'param 2 is not of type int or long'

  npar1 = n_elements(param1)
  npar2 = n_elements(param2)

  IF npar1 NE npar2 THEN message,'# in param1 ne # in param2'
  indices = lonarr(npar1)

  hist = histogram(param1, rever=rev_ind)

  whist = where(hist NE 0, nhist)
  
  beg = 0L
  FOR i=0L, nhist-1 DO BEGIN 

      binnum = whist[i]

      w=rev_ind( rev_ind[binnum]:rev_ind[binnum+1]-1 )
      nw = n_elements(w)

      sw=sort(param2[w])

      indices[beg:beg+nw-1] = w[sw]

      beg = beg+nw

  ENDFOR 
      
  return,indices

END 
