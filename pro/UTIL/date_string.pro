FUNCTION date_string, dtstring

  IF keyword_set(help) THEN BEGIN
      print,'-Syntax: result = date_string(help=help)'
      return,''
  ENDIF 

  CALDAT, systime(/julian), Month, Day, Year, Hour, Minute, Second
  dtstring = ntostr(day)+'-'+ntostr(month)+'-'+ntostr(year)

  return,dtstring
END
