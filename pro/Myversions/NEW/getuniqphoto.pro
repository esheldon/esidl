FUNCTION getuniqphoto, run, rerun, camcol, field, id, nuniq

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    GETUNIQPHOTO
;       
; PURPOSE:
;    Find uniq objects by run,rerun,camcol,field,id, which uniquely
;    identify PHOTO output objects
;
; CALLING SEQUENCE:
;    uniq = getuniqphoto(run, rerun, camcol, field, id [,num])
;
; INPUTS: 
;    run, rerun, camcol, field, id
;
; OPTIONAL INPUTS:
;    none
;
; KEYWORD PARAMETERS:
;    none
;       
; OUTPUTS: 
;    uniq: the unique elements.
;
; OPTIONAL OUTPUTS:
;    num: number of unique elements
;
; CALLED ROUTINES:
;    REM_DUP
; 
; PROCEDURE: 
;    make a long64 word from the input integer values and find the
;    unique elements
;	
;
; REVISION HISTORY:
;    Created: 21-9-2000  Erin Scott Sheldon  UofM
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() EQ 0 THEN BEGIN 
     print,'-Syntax: uniq = getuniqphoto(run, rerun, camcol, field, id [,num])'
     print,''
     print,'Use doc_library,"getuniqphoto"  for more help.'  
     return,-1
  ENDIF 

  ten=ulong64(10)

  super=ulong64(id)+ulong64(field)*ten^6 +ulong64(camcol)*ten^11 $
    +ulong64(rerun)*ten^13 + ulong64(run)*ten^16

  uniq = rem_dup(super)
  nuniq = n_elements(uniq)

  return,uniq
END 
