FUNCTION geterrtag, struct

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: result = geterrtag(struct)'
      return,-1
  ENDIF


  IF NOT tag_exist(struct,'momerr',index=momerr) THEN BEGIN
      IF NOT tag_exist(struct,'uncert',index=momerr) THEN BEGIN 
          message,'Struct must have "momerr" or "uncert" tag'
      ENDIF 
  ENDIF 

  return, momerr

END 
