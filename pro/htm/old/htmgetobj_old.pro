
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    htmGetObj
;       
; PURPOSE:
;    Get the objects from an htm tree leaf.
;    Currently objects consist of a simple structure containing
;    only an index and a pointer to the next object. This returns the index.
;    of all objects.
; 
;    HTM = hierarchical triangular mesh see http://www.sdss.jhu.edu/
;
; CALLING SEQUENCE:
;    objects = htmGetObj(id, tree, count)
;
; INPUTS: 
;    id: the leaf node id.
;    tree: the HTM tree, as build by htmBuildTree
;
; OPTIONAL INPUTS:
;    NONE
;
; KEYWORD PARAMETERS:
;    NONE
;       
; OUTPUTS: 
;    The objects (only indices for now)
;
; OPTIONAL OUTPUTS:
;    count: the number of objects found
;
; CALLED ROUTINES:
;    htmGetIdNode
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


FUNCTION htmGetObj_old, id, tree, count

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: objects = htmGetObj(id, tree [, count])'
      print,''
      print,'Use doc_library,"htmGetObj"  for more help.'  
      return,-1
  ENDIF 

  idnode = htmGetIdNode(id, tree)

  count = (*idnode).nobj
  IF count EQ 0 THEN return,-1

  indices = replicate(ulong(0), count)
  tmp = (*idnode).obj
  indices[0] = (*tmp).index
  next = (*tmp).nextobj
  FOR i=1L, count-1 DO BEGIN 
      tmp = next
      indices[i] = (*tmp).index
      next = (*tmp).nextobj
  ENDFOR 
  return,indices

END 
