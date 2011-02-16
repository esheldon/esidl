
PRO db_struct2def, struct, names, typedefs

  IF n_elements(struct) NE 1 THEN message,'Just send a single struct'

  tags = strlowcase( tag_names(struct) )
  ntags = n_elements(tags)

  names = strarr(ntags)
  typedefs = names
  comments = names

  
  FOR i=0L, ntags-1 DO BEGIN 
      t = size(struct.(i), /struct)
      
      names[i] = tags[i]

      IF (t.n_dimensions) GT 0 THEN BEGIN 
          ;; array
          dim=t.dimensions[0:t.n_dimensions-1]
          addst = '('+strjoin(ntostr(dim),',')+')'
          names[i] = names[i] + addst
      ENDIF 

      typedefs[i] = db_idl2type(t)

  ENDFOR 

END 
