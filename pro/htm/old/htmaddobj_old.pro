
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    htmAddObj
;       
; PURPOSE:
;    Add objects to a HTM tree leaf (hierarchical triangular mesh).
;    Currently objects consist of a simple structure containing
;    only an index and a pointer to the next object.
;  
;    HTM = hierarchical triangular mesh see http://www.sdss.jhu.edu/
;
; CALLING SEQUENCE:
;    htmAddObj, id, indices, tree
;
; INPUTS: 
;    id: the leaf node id.
;    indices: the indices of objects. e.g. indices in a structure
;             containing ra,dec
;    tree: the HTM tree, as build by htmBuildTree
;
; OPTIONAL INPUTS:
;    NONE
;
; KEYWORD PARAMETERS:
;    NONE
;       
; OUTPUTS: 
;    the augmented tree.
;
; OPTIONAL OUTPUTS:
;    NONE
;
; COMMON BLOCKS:
;     COMMON htmcomm, node, object, base, true, false
;
; CALLED ROUTINES:
;    htmGetIdNode
;    htmGetLastObj
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    11-DEC-2000 Creation Erin Scott Sheldon UofMich
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO htmAddObj_old, id, indices, tree

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: htmAddObj, ids, indices, tree'
      print
      print,'Use doc_library,"htmAddObj" for more help.'
      return
  ENDIF 

  COMMON htmcomm, node, object, base, true, false
  IF n_elements(node) EQ 0 THEN htmsetup

  ; id is the node id
  ; indices are labels, such as the indices of objects
  ; in an idl structure.

  nind = n_elements(indices)
  starti = 0L

  idnode = htmGetIdNode(id, tree)
  nobj = (*idnode).nobj
  IF nobj NE 0 THEN BEGIN
      print,'Finding Last Object'
      current = htmGetLastObj( (*idnode).obj, nobj )
  ENDIF ELSE BEGIN  
      print,'No objects. Starting New Linked List'
      (*idnode).obj = ptr_new(object)
      current = (*idnode).obj    
      (*current).index = indices[0]
      starti = starti + 1L
  ENDELSE 
  (*idnode).nobj = nobj+nind

  next = ptr_new(object)
  FOR i=starti, nind-1 DO BEGIN 
      (*current).nextobj = next

      tmp = (*current).nextobj
      (*tmp).index = indices[i]
      IF i NE nind-1 THEN BEGIN 
          next = ptr_new(object)
          current = tmp
      ENDIF 
  ENDFOR 


END 
