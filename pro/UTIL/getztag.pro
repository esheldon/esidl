FUNCTION getztag, struct, zerrtag=zerrtag

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: result = getztag(struct)'
      return,-1
  ENDIF

  IF NOT tag_exist(struct[0],'z1d',index=ztag) THEN BEGIN
      IF NOT tag_exist(struct[0],'z',index=ztag) THEN BEGIN 
          IF NOT tag_exist(struct[0],'photoz_z',index=ztag) THEN BEGIN
              ;;print,'Lens struct does not have "Z" or "PHOTOZ" or "Z1D" flag'
          ENDIF 
      ENDIF 
  ENDIF 

  IF NOT tag_exist(struct[0],'z1d_error',index=zerrtag) THEN BEGIN 
      IF NOT tag_exist(struct[0],'zerr',index=zerrtag) THEN BEGIN 
          IF NOT tag_exist(struct[0],'photoz_zerr',index=zerrtag) THEN BEGIN 
              ;;print,'Lens struct does not have "Z" or "PHOTOZ" or "Z1D" flag'
          ENDIF 
      ENDIF 
  ENDIF

  return, ztag

END
