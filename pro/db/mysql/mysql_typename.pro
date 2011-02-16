FUNCTION mysql_typename, mysql_type

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: typename = mysql_typename(mysql_type)'
      return,''
  ENDIF 

  ;; These are the names used for declaration of a table

  CASE mysql_type OF
      0: return,'DECIMAL'
      1: return,'TINYINT'
      2: return,'SMALLINT'
      3: return,'INT'
      4: return,'FLOAT'
      5: return,'DOUBLE'
      6: return,'NULL'
      7: return,'TIMESTAMP'
      8: return,'BIGINT'
      9: return,'MEDIUMINT'
      10: return,'DATE'
      11: return,'TIME'
      12: return,'DATETIME'
      13: return,'YEAR'
      14: return,'NEWDATE'
      247: return,'ENUM'
      248: return,'SET'
      249: return,'TINY_BLOB'
      250: return,'MEDIUM_BLOB'
      251: return,'LONG_BLOB'
      252: return,'BLOB'
      253: return,'VAR_STRING'
      254: return,'STRING'
      ELSE: BEGIN 
          message,'Unknown mysql type number: '+ntostr(mysql_type),/inf
          return,'UNKNOWN'
      END 
  ENDCASE 

  ;; These are the defs in the .h file. I decided to use the above for
  ;; clarity 
  CASE mysql_type OF
      0: return,'DECIMAL'
      1: return,'TINY'
      2: return,'SHORT'
      3: return,'LONG'
      4: return,'FLOAT'
      5: return,'DOUBLE'
      6: return,'NULL'
      7: return,'TIMESTAMP'
      8: return,'LONGLONG'
      9: return,'INT24'
      10: return,'DATE'
      11: return,'TIME'
      12: return,'DATETIME'
      13: return,'YEAR'
      14: return,'NEWDATE'
      247: return,'ENUM'
      248: return,'SET'
      249: return,'TINY_BLOB'
      250: return,'MEDIUM_BLOB'
      251: return,'LONG_BLOB'
      252: return,'BLOB'
      253: return,'VAR_STRING'
      254: return,'STRING'
      ELSE: BEGIN 
          message,'Unknown mysql type number: '+ntostr(mysql_type),/inf
          return,'UNKNOWN'
      END 
  ENDCASE 
END 
