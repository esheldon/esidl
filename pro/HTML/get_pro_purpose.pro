PRO get_pro_purpose, file, purp_line

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: get_pro_purpose, file, purp_line'
      return
  ENDIF 

  purp_line = ''
  openr, lun, file, /get_lun

  line = ''
  cont=1
  WHILE (NOT EOF(lun)) AND cont DO BEGIN 

      readf, lun, line
      
      IF strpos(line,'PURPOSE') NE -1 THEN BEGIN 
          ;; read next line for hint at purpose
          readf, lun, purp_line
          ;; remove beginning semicolon and space
          len = strlen(purp_line)
          IF strmid(purp_line, 0, 1) EQ ';' THEN BEGIN 
              purp_line = strmid(purp_line, 1,len-1)
              len = strlen(purp_line)
          ENDIF 
          purp_line = strtrim(string(purp_line), 2)
          cont=0
      ENDIF 

  ENDWHILE 

  free_lun, lun

END 
