PRO invbw

  tvlct,r,g,b,/get
  IF r[0] NE 0 THEN BEGIN 
      loadct,0
      !p.color=255
  ENDIF ELSE BEGIN 
      loadct,41,file="$IDL_DIR/resource/colors/colors1.tbl"
      !P.COLOR=255
  ENDELSE 
  
END 
