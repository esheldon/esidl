PRO ppercent, percent

  back=string(8B)
  pstr = ntostr( long(percent) )
  IF strlen(pstr) EQ 1 THEN pstr = '0'+pstr

  print,back+back+back+pstr+'%',format='(a,$)'

END 
