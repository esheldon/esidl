FUNCTION db_idl2type, szStruct, string_length=string_length

  CASE szStruct.type_name OF
      'BYTE':    return, 'B*1'
      'INT':     return, 'I*2'
      'UINT':    return, 'U*2'
      'LONG':    return, 'I*4'
      'ULONG':   return, 'U*4'
      'LONG64':  return, 'I*8'
      'ULONG64': return, 'U*8'

      'FLOAT':   return, 'R*4'
      'DOUBLE':  return, 'R*8'

      'STRING':  BEGIN 
          IF n_elements(string_length) EQ 0 THEN BEGIN 
              message,'string_length not entered, using 20',/inf
              string_length=20
          ENDIF 
          return, 'C*'+ntostr(string_length)
      END 
  ENDCASE 

END 
