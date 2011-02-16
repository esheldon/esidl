PRO pred, string, bold=bold

  redout='31'
  IF keyword_set(bold) THEN BEGIN 
      spawn,'echo "\033['+redout+';01m'+string+'\033[0;00m"'
  ENDIF ELSE BEGIN 
      spawn,'echo "\033['+redout+';'+redout+'m'+string+'\033[0;00m"'
  ENDELSE 

END 
