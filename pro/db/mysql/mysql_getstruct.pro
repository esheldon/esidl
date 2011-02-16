FUNCTION mysql_getstruct, field_names, $
                          mysql_field_types, $
                          mysql_field_flags, $
                          nvalues=nvalues

  ;; Get the data types from the field definitions
  ;; and create a structure

  nf  = n_elements(field_names)
  nft = n_elements(mysql_field_types)
  nff = n_elements(mysql_field_flags)
  nv  = n_elements(nvalues)

  IF nft NE nf THEN BEGIN 
      message,'field types must be same size as field names',/inf
      return,-1
  ENDIF 
  IF nff NE nf THEN BEGIN 
      message,'field flags must be same size as field names',/inf
      return,-1
  ENDIF 
  
  IF nv NE nf THEN BEGIN 
      IF nv EQ 0 THEN BEGIN 
          nvalues = replicate(1, nf)
      ENDIF ELSE BEGIN 
          message,'nvalues must be same size as field names',/inf
          return,-1
      ENDELSE 
  ENDIF 

  idl_field_types = strarr(nf)

  FOR i=0L, nf-1 DO BEGIN 
      idl_field_types[i] = $
        mysql_idltype(mysql_field_types[i],mysql_field_flags[i], $
                      nvalues=nvalues[i])
  ENDFOR 

  struct = mrd_struct(field_names, idl_field_types, 1)
  return,struct

END 
