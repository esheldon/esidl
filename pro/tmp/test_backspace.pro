
PRO test_backspace, w


  sequence = ['-','\','|', '/']

  back = string(8B)

  IF w EQ 1 THEN BEGIN 
      j=0
      WHILE 1 DO BEGIN
          IF j EQ 4 THEN j=0
          str = ntostr(j)
          print,back,format='(a,$)'
          print,sequence[j],format='(a,$)'
          j=j+1
          wait,1
      ENDWHILE 
  ENDIF 

  IF w EQ 2 THEN BEGIN 
      
      ntot = 100
      ntotstr = field2string(ntot)
      
      nstr=9
      back=string(8b)
      FOR i=1,nstr-1 DO back=back+back
      lasti = '0000/'
      str = lasti+ntotstr
      FOR i=0L, ntot-1 DO BEGIN 
          
          istr = field2string(i+1)
          str = repstr(str, lasti, istr+'/')
          lasti = istr+'/'

          print,back,format='(a,$)'
          print,str,format='(a,$)'
          wait,1
      ENDFOR 
  ENDIF 
  print

  return
END 
