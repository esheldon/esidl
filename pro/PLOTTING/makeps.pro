PRO makeps, name, noland=noland

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: makeps
;       
; PURPOSE: wrapper for begplot
;	
;
; CALLING SEQUENCE: makeps, name, noland=noland
;
; INPUTS: name:  name of ps file
;
; INPUT KEYWORD PARAMETERS: noland:  if set, not landscape
;
;
; REVISION HISTORY:
;	Erin Scott Sheldon  UMich  6/17/99
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if N_params() eq 0 then begin
      print,'-Syntax: makeps, name, noland=noland'
      print,''
      print,'Use doc_library,"makeps"  for more help.'  
      return
  endif
  
  IF keyword_set(noland) THEN landscape = 0 ELSE landscape = 1
  begplot,name=name,/invbw,landscape=landscape


return
end































