pro esrextract, param_struct, filename
;+
; NAME:
;       REXTRACT
; PURPOSE:
;	Run sextractor on an image using the parameters passed in param_struct
;
; CALLING SEQUENCE:
;       rextract,param_struct, filename
;
; INPUTS:
;	param_struct: a structure containing all the parameters needed
;		      to run sextractor
;	filename: image to be processed
;       
; OUTPUTS:
;	
;
; OPTIONAL OUTPUT ARRAYS:
;
; INPUT KEYWORD PARAMETERS:
; 
; PROCEDURE: This processes an image on disk using the parameters provided.
;	
;
; REVISION HISTORY:
;	Tim McKay	UM	1/8/98
;-
 On_error,2              ;Return to caller

 if N_params() ne 2 then begin
        print,'Syntax - esrextract, param_struct, filename
        return
 endif

 tags=tag_names(param_struct) 
 tagsize=size(tags)
 
 if tagsize(1) ne 35 then begin
	print, 'These are not valid parameters!!!'
	return
 endif
 
 cmd_string = '/sdss/products/sextractor/sextractor1.2b10b/source/sex '
 cmd_string = cmd_string+filename
 cmd_string = cmd_string+' -c /sdss/products/sextractor/sextractor1.2b10b/config/rotse.sex '
 

 for i = 0,33,1 do begin
  	tag = tags(i)
	val=string(param_struct.(i))
	cmd_string = cmd_string+'-'+tag+' '+val+' '
 endfor
   
 print, cmd_string
 spawn, cmd_string

 return

 end








