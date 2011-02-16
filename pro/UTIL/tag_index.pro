FUNCTION tag_index, struct, tag_name

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: tagid = tag_index(struct, tag_name)'
      return,-1
  ENDIF 

  tmp = tag_exist(struct, tag_name, index=index)
  return,index

END 
