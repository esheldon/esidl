
function mrd_struct2, names, values, nrow, $
    structyp=structyp, tempdir=tempdir,silent=silent,old_struct=old_struct
;+
; NAME:
;       MRD_STRUCT
; PURPOSE:
;       Return a structure as defined in the names and values data.
;        !! The idlastron guys put this new stuff in their code, so
;        there is no need for this procedure.
; CALLING SEQUENCE:
;       struct = MRD_STRUCT(NAMES, VALUES, NROW,                $
;                   STRUCTYP=structyp,                          $
;                   TEMPDIR=tempdir)
; INPUT PARAMETERS:
;       NAMES   = A string array of names of structure fields.
;       VALUES  = A string array giving the values of the structure
;                 fields.  See examples below.
;       NROW    = The number of elements in the structure array.
;       
; RETURNS:
;       A structure as described in the parameters or 0 if an error
;       is detected.
;
; OPTIONAL KEYWORD PARAMETERS:
;       STRUCTYP = The structure type.  Since IDL does not allow the
;                  redefinition of a named structure it is an error
;                  to call MRD_STRUCT with different parameters but
;                  the same STRUCTYP in the same session.  If this
;                  keyword is not set an anonymous structure is created.
;       TEMPDIR  = Just left in so won't conflict with calls in mrdfits
;                       E.S.S.
; COMMON BLOCKS:
;	None
;
; PROCEDURE:
;       A structure definition is created using the parameter values.
; EXAMPLES:
;       str = mrd_struct(['fld1', 'fld2'], ['0','dblarr(10,10)'],3)
;       print, str(0).fld2(3,3)
;
;       str = mrd_struct(['a','b','c','d'],['1', '1.', '1.d0', "'1'"],1)
;               ; returns a structure with integer, float, double and string
;               ; fields.
; MODIFICATION HISTORY:
;       Created by T. McGlynn October, 1994.
;       Modified by T. McGlynn September, 1995.
;          Added capability to create substructures so that structure
;          may contain up to 4096 distinct elements.  [This can be
;          increased by futher iteration of the process used if needed.]
;	Converted to IDL V5.0   W. Landsman   September 1997
;	Removed V4.0 reference to common block  October 1997
;       Completely rewrote Erin Scott Sheldon 28-OCT-2000 U. of Michigan
;-


; Check that the number of names is the same as the number of values.

  nel = n_elements(names)
  IF nel NE n_elements(values) THEN return, 0

  ;;;;;;;;;;;;;;;;;;;;;;;;
  ;; structure name
  ;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(structyp) EQ 0 THEN struct_type = '' ELSE struct_type=structyp

  ;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; create first element
  ;;;;;;;;;;;;;;;;;;;;;;;;;

  command = 'strtmp = create_struct("'+names[0]+'", '+values[0]+')'
  res=execute(command)
  IF res EQ 0 THEN return,0     ;Check for errors

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Add any more elements required
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR i=1L, nel-1 DO BEGIN
      command = 'strtmp = create_struct(strtmp, "'+names[i]+'", '+values[i]+')'
      res = execute(command)
      IF res EQ 0 THEN return,0     ;Check for errors
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; see if we need to return an array of structs
  ;; also check if needs to be a named structure
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF nrow GT 1 THEN BEGIN 
      command = 'struct = replicate(create_struct(name="'+$
        struct_type+'", strtmp), '+strtrim(long(nrow),2)+ ')'
  ENDIF ELSE BEGIN 
      IF struct_type EQ '' THEN return,strtmp $
      ELSE command = 'struct = create_struct(name="'+struct_type+'", strtmp)'
  ENDELSE 
  res = execute(command)
  IF res EQ 0 THEN return,0 ELSE return,struct     ;Check for errors/return

END 
