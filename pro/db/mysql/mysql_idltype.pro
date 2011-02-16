FUNCTION mysql_idltype, mysql_type, mysql_flags, nvalues=nvalues

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: type = mysql_idltype(mysql_type, mysql_flags, nvalues=)'
      return,''
  ENDIF 

  FIELD_IS_UNSIGNED_FLAG = 32

  is_unsigned = ( (mysql_flags AND FIELD_IS_UNSIGNED_FLAG) NE 0)
  
  ;; Allow the type to be an array
  IF n_elements(nvalues) EQ 0 THEN nvalues = 1

  IF nvalues[0] EQ 1 THEN BEGIN 

      CASE mysql_type OF

          ;; tinyint No idl type like this when signed, so use int
          1: IF is_unsigned THEN tagtype = '0b' else tagtype='0'
          ;; int
          2: IF is_unsigned THEN tagtype = '0u' ELSE tagtype = '0'
          ;; long
          3: IF is_unsigned THEN tagtype='0UL' else tagtype = '0L'
          ;; float
          4: tagtype = '0.0'
          ;; double
          5: tagtype = '0d'
          ;; long long
          8: IF is_unsigned THEN tagtype='0ULL' else tagtype = '0LL'
          ;; int 24, just use long
          9: IF is_unsigned THEN tagtype='0UL' else tagtype = '0L'
          ;; Date
          10: tagtype='""'
          ;; Time
          11: tagtype = '""'
          ;; Date Time
          12: tagtype = '""'
          ;; Year
          13: tagtype = '""'
          ;; Newdate
          14: tagtype = '""'
          ;; Enum: be safe and use long
          247: tagtype = '0L'
          ;; Var String
          253: tagtype='""'
          ;; String
          254: tagtype='""'
          ELSE: BEGIN 
              typename = mysql_typename(mysql_type)
              message,'Unsupported mysql type: '+typename
          END 
      ENDCASE 

  ENDIF ELSE BEGIN 

      CASE mysql_type OF

          ;; tinyint No idl type like this when signed, so use int
          1: IF is_unsigned THEN tagtype = 'bytarr' else tagtype='intarr'
          ;; int
          2: IF is_unsigned THEN tagtype = 'uintarr' ELSE tagtype = 'intarr'
          ;; long
          3: IF is_unsigned THEN tagtype='ulonarr' else tagtype = 'lonarr'
          ;; float
          4: tagtype = 'fltarr'
          ;; double
          5: tagtype = 'dblarr'
          ;; long long
          8: IF is_unsigned THEN tagtype='ulon64arr' else tagtype = 'lon64arr'
          ;; int 24, just use long
          9: IF is_unsigned THEN tagtype='ulonarr' else tagtype = 'lonarr'
          ;; Date
          10: tagtype='strarr'
          ;; Time
          11: tagtype = 'strarr'
          ;; Date Time
          12: tagtype = 'strarr'
          ;; Year
          13: tagtype = 'strarr'
          ;; Newdate
          14: tagtype = 'strarr'
          ;; Enum: be safe and use long
          247: tagtype = 'lonarr'
          ;; Var String
          253: tagtype='strarr'
          ;; String
          254: tagtype='strarr'
          ELSE: BEGIN 
              typename = mysql_typename(mysql_type)
              message,'Unsupported mysql type: '+typename
          END 
      ENDCASE 
      nstr = '(' + ntostr(nvalues) +')'
      tagtype = tagtype + nstr

  ENDELSE 


  return,tagtype

END 
