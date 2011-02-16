PRO getstructlambda, struct, clambda, ceta, radec=radec

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: getstructlambda, struct, clambda, ceta, /radec'
      print,' if /radec it finds ra,dec tags: getstructlambda,struct,ra,dec,/radec'
      return
  ENDIF 

  IF tag_exist(struct[0], 'ra') AND tag_exist(struct[0], 'dec') THEN BEGIN
      ra_exist=1 
  ENDIF ELSE ra_exist=0

  IF tag_exist(struct[0],'clambda') AND tag_exist(struct[0],'ceta') THEN BEGIN
      clambda_exist=1 
  ENDIF ELSE clambda_exist=0

  IF NOT keyword_set(radec) THEN BEGIN 
      IF NOT clambda_exist THEN BEGIN 
          IF NOT ra_exist THEN $
            message,'Neither (CLAMBDA, CETA) or (RA,DEC) found in lens tags'
          eq2csurvey, struct.ra, struct.dec, clambda, ceta
      ENDIF ELSE BEGIN 
          clambda = struct.clambda
          ceta = struct.ceta
      ENDELSE 
  ENDIF ELSE BEGIN 
      IF NOT ra_exist THEN BEGIN 
          IF NOT clambda_exist THEN $
            message,'Neither (CLAMBDA, CETA) or (RA,DEC) found in lens tags'
          csurvey2eq, struct.clambda, struct.ceta, clambda, ceta
      ENDIF ELSE BEGIN 
          clambda = struct.ra
          ceta = struct.dec
      ENDELSE 

  ENDELSE 

END 
