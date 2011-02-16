FUNCTION htmAddNode_old, lowerid, addlevel, addbit

  COMMON htmcomm, node, object, base, true, false
  
  tmp = ptr_new(node)
  IF addbit THEN id = lowerid + base^addlevel ELSE id = lowerid
  (*tmp).id = id

  IF addlevel NE 0 THEN BEGIN 
      nextlevel = addlevel-1
      sendlevel = nextlevel-1
      sub1id = id
      (*tmp).sub1 = htmAddNode_old(sub1id, sendlevel, false)
      sub2id = id
      (*tmp).sub2 = htmAddNode_old(sub2id, sendlevel, true)
      sub3id = id + base^nextlevel
      (*tmp).sub3 = htmAddNode_old(sub3id, sendlevel, false)
      sub4id = sub3id
      (*tmp).sub4 = htmAddNode_old(sub4id, sendlevel, true)
  ENDIF ELSE print,id
  return,tmp

END 
