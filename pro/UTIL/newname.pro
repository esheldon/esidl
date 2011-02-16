FUNCTION newname, testname

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: newn = newname(oldname)'
      print
      print,'Names should be in this form: stuff_N1.ext, or N2,N3,etc'
      print,'newname will increment the number.'
      return,""
  ENDIF 

  result = str_sep(testname, '.')
  
  str = '_N'
  a = rstrpos(result[0], str)
  
  IF a EQ -1 THEN BEGIN
      print,'Improper string: ',testname
      return,testname
  ENDIF ELSE BEGIN
      beg = strmid(result[0],0,a)
      nr = strmid(result[0],a+2)
      eend = '_N'+ntostr( long(nr)+1 )
      
      newname = beg + eend + '.'+result[1]
      return,newname
  ENDELSE 
    

END 
