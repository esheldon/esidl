PRO htmFill2Leaf_old, bnode, level, leafid

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: htmFill2Leaf, node, level, leafid'
      return
  ENDIF 

  COMMON htmcomm, node, object, base, true, false

  id = (*bnode).id
  tlevel = level+1
  nextlevel = tlevel-1
  sendlevel = nextlevel-1

  ;; figure out which way to go
  tsub1id = id
  tsub2id = id + base^sendlevel
  tsub3id = id + base^nextlevel
  tsub4id = tsub3id + base^sendlevel
  
  IF tlevel NE 0 THEN BEGIN 
      CASE 1 OF 
          (leafid GE tsub4id): BEGIN 
              sub4id = id + base^nextlevel
              (*bnode).sub4 = htmAddNode(leafid, sub4id, sendlevel, true)
          END
          (leafid GE tsub3id): BEGIN 
              sub3id = id + base^nextlevel
              (*bnode).sub3 = htmAddNode(leafid, sub3id, sendlevel, false)
          END 
          (leafid GE tsub2id): BEGIN
              sub2id = id
              (*bnode).sub2 = htmAddNode(leafid, sub2id, sendlevel, true)
          END 
          ELSE: BEGIN 
              sub1id = id
              (*bnode).sub1 = htmAddNode(leafid, sub1id, sendlevel, false)
          END 
      ENDCASE 
  ENDIF 

END 
