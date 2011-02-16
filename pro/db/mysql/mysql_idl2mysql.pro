FUNCTION mysql_idl2mysql, tname, length=length
  
  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: mysql_typename = mysql_idl2mysql(type_name, length=)'
      print,'You must enter a length for string types'
      return,''
  ENDIF 

  ;; Convert an idl type name to a mysql type name
  CASE strupcase(tname) OF 

      'BYTE': return,'TINYINT UNSIGNED'
      'INT': return,'SMALLINT'
      'UINT': return,'SMALLINT UNSIGNED'
      'LONG': return,'INT'
      'ULONG': return,'INT UNSIGNED'
      'LONG64': return,'BIGINT'
      'ULONG64': return,'BIGINT UNSIGNED'
      'FLOAT': return,'FLOAT'
      'DOUBLE': return,'DOUBLE'
      'STRING': BEGIN 
          IF n_elements(length) EQ 0 THEN BEGIN 
              message,'You must enter a length for string types'
          ENDIF 
          return,'VARCHAR('+strtrim( long(length[0]) ,2)+')'
      END 
  ENDCASE 

END 
