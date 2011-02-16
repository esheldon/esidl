FUNCTION corrselect_flag, type

  IF size(type, /tname) EQ 'STRING' THEN BEGIN 
      CASE strupcase(type) OF
          'GOODU': return, 2b^0
          'GOODG': return, 2b^1
          'GOODR': return, 2b^2
          'GOODI': return, 2b^3
          'GOODZ': return, 2b^4
          ELSE: message,'Unrecognized flag: '+strupcase(type)
      ENDCASE
  ENDIF ELSE BEGIN 
      CASE fix(type) OF
          0: return, 2b^0
          1: return, 2b^1
          2: return, 2b^2
          3: return, 2b^3
          4: return, 2b^4
          ELSE: message,'Unrecognized flag: '+ntostr(type)
      ENDCASE
  ENDELSE 
END 
