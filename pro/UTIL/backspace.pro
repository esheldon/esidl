FUNCTION backspace, num
  back = string(8B)
  IF n_elements(num) NE 0 THEN BEGIN 
      back = mkstr(num, val=back)
  ENDIF 
  return, back
END 
