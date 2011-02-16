
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    htmBuildTree
;       
; PURPOSE:
;    Build an HTM tree to the specified depth.
;    HTM = hierarchical triangular mesh. see http://www.sdss.jhu.edu/
;
; CALLING SEQUENCE:
;    htmBuildTree, depth, tree
;
; INPUTS: 
;    depth: The depth of the tree. 0th level splits the sphere into 8 
;           triangles. Each triangle is split into 4 at level 1, and
;           each of these is split into 4 at the next level, etc.
;           level 8 is the practical limit in idl.
;
; OPTIONAL INPUTS:
;    NONE
;
; KEYWORD PARAMETERS:
;    NONE
;       
; OUTPUTS: 
;    tree: an HTM tree to specified depth.
;
; OPTIONAL OUTPUTS:
;    NONE
;
; COMMON BLOCKS:
;     COMMON htmcomm, node, object, base, true, false
;
; CALLED ROUTINES:
;    htmAddNode (recursive)
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

PRO htmBuildTree_old, depth, tree

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: htmBuildTree, depth, tree'
      print,''
      print,'Use doc_library,"htmBuildTree"  for more help.'  
      return
  ENDIF 

  COMMON htmcomm, node, object, base, true, false

  IF n_elements(node) EQ 0 THEN htmsetup

  base = ulong(2)
  true = 1b
  false = 0b

  ;; set up bottom level
  bottom_level = htmBottomLevel(depth)
  addlevel = bottom_level-1
  bottomid = base^bottom_level
  maxid = base^(bottom_level+1)-1
  print,'bottom level',bottom_level
  print,'bottom id = ',bottomid
  print,'addlevel = ',addlevel
  print,'maxid = ',maxid

  print
  tree = ptr_new(node)
;  (*tree).id = bottomid

  ;; north
  (*tree).sub4 = htmAddNode_old(bottomid, addlevel, true )
  ;; south
  (*tree).sub3 = htmAddNode_old(bottomid, addlevel, false)

END 
