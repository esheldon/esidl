FUNCTION htmGetIdNode_old, id, tree

  IF n_params() LT 2 THEN BEGIN
      message,'incorrect # of args',/cont
      message,'-Syntax: node=htmGetIdNode_old(id, tree)'
  ENDIF 

  COMMON htmcomm, node, object, base, true, false
  IF n_elements(node) EQ 0 THEN htmsetup

  tmp = (*tree).sub4
  IF ptr_valid(tmp) THEN BEGIN 
      IF id GE (*tmp).id THEN return,htmGetIdNode_old(id,tmp)
  ENDIF 
  tmp = (*tree).sub3
  IF ptr_valid(tmp) THEN BEGIN 
      IF id GE (*tmp).id THEN return,htmGetIdNode_old(id,tmp)
  ENDIF 
  tmp = (*tree).sub2
  IF ptr_valid(tmp) THEN BEGIN 
      IF id GE (*tmp).id THEN return,htmGetIdNode_old(id,tmp)
  ENDIF 
  tmp = (*tree).sub1
  IF ptr_valid(tmp) THEN BEGIN 
      IF id GE (*tmp).id THEN return,htmGetIdNode_old(id,tmp)
  ENDIF 
  ;; if ptr not valid for any, then we are done
  return,tree

END 
