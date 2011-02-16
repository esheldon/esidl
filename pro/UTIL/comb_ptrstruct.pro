PRO comb_ptrstruct, base_struct, ptrlist, numlist, outstruct,nofree=nofree

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: comb_ptrstruct, base_struct, ptrlist, numlist, outstruct, /nofree'
      print
      print,'Pointer list is freed after copying into structure unless /nofree'
      return
  ENDIF 

  nf = n_elements(numlist)
  IF n_elements(ptrlist) NE nf THEN message,'ptrlist must be same size as numlist'

  ;; Can't use total() because it converts to double
  ntot = 0L
  FOR i=0L, nf-1 DO ntot = ntot + numlist[i]
  outstruct = replicate(base_struct, ntot )

  beg=0L
  FOR fi=0L, nf-1 DO BEGIN 

      IF numlist[fi] NE 0 THEN BEGIN 
          outstruct[beg:beg+numlist[fi]-1] = *ptrlist[fi]
      ENDIF 
      IF NOT keyword_set(nofree) THEN ptr_free, ptrlist[fi]
      beg=beg+numlist[fi]
  
  ENDFOR 

  return
END 
