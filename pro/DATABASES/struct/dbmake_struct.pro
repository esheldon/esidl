FUNCTION dbmake_struct, items, idltype, numvals, nrows=nrows

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: struct = dbmake_struct(items, idltype, numvals, nrows=)'
      return,-1
  ENDIF 

  nitem = n_elements(items)

  values = strarr(nitem)

  FOR i=0L, nitem-1 DO BEGIN 

      CASE idltype[i] OF
          1: BEGIN 
              IF numvals[i] EQ 1 THEN values[i] = '0b' $
              ELSE  values[i] = 'BYTARR('+ntostr(numvals[i])+')'
          END 
          2: BEGIN 
              IF numvals[i] EQ 1 THEN values[i] = '0' $
              ELSE  values[i] = 'INTARR('+ntostr(numvals[i])+')'
          END 
          3: BEGIN 
              IF numvals[i] EQ 1 THEN values[i] = '0L' $
              ELSE  values[i] = 'LONARR('+ntostr(numvals[i])+')'
          END 
          4: BEGIN 
              IF numvals[i] EQ 1 THEN values[i] = '0.0' $
              ELSE  values[i] = 'FLTARR('+ntostr(numvals[i])+')'
          END 
          5: BEGIN 
              IF numvals[i] EQ 1 THEN values[i] = '0D' $
              ELSE  values[i] = 'DBLARR('+ntostr(numvals[i])+')'
          END 
          7: BEGIN 
              values[i] = strjoin(strarr(numvals[i]))
          END 
          12: BEGIN 
              IF numvals[i] EQ 1 THEN values[i] = '0U' $
              ELSE  values[i] = 'UINTARR('+ntostr(numvals[i])+')'
          END 
          13: BEGIN 
              IF numvals[i] EQ 1 THEN values[i] = '0UL' $
              ELSE  values[i] = 'ULONARR('+ntostr(numvals[i])+')'
          END 
          14: BEGIN 
              IF numvals[i] EQ 1 THEN values[i] = '0LL' $
              ELSE  values[i] = 'LON64ARR('+ntostr(numvals[i])+')'
          END 
          15: BEGIN 
              IF numvals[i] EQ 1 THEN values[i] = '0ULL' $
              ELSE  values[i] = 'ULON64ARR('+ntostr(numvals[i])+')'
          END 
          ELSE: message,'Unsupported type : ',idltype[i]

      ENDCASE 
  ENDFOR 
  
  IF n_elements(nrows) EQ 0 THEN nrows=1
  struct = mrd_struct(items, values, nrows)
  return,struct
END 
  
