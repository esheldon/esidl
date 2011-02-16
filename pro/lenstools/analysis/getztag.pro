FUNCTION getztag, struct, zerrtag=zerrtag

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: result = getztag(struct)'
      return,-1
  ENDIF

  IF NOT tag_exist(struct,'z1d',index=ztag) THEN BEGIN
      IF NOT tag_exist(struct,'z',index=ztag) THEN BEGIN 
          IF NOT tag_exist(struct,'photoz_z',index=ztag) THEN BEGIN
              message,'Lens structure must have "Z" or "PHOTOZ" or "Z1D" flag',/inf
              return,-1
          ENDIF 
      ENDIF 
  ENDIF 

  IF NOT tag_exist(struct,'z1d_error',index=zerrtag) THEN BEGIN 
      IF NOT tag_exist(struct,'zerr',index=zerrtag) THEN BEGIN 
          IF NOT tag_exist(struct,'photoz_zerr',index=zerrtag) THEN BEGIN 
          ENDIF 
      ENDIF 
  ENDIF

  return, ztag

END
