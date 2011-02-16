pro call_sex, param_struct, filename
;+
; NAME:
;       REXTRACT
; PURPOSE:
;	Run sextractor on an image using the parameters passed in param_struct
;
; CALLING SEQUENCE:
;       deep_extract,param_struct, filename
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
        print,'Syntax - call_sex, param_struct, filename'
        return
 endif

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; see which tags are defined
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 tags=tag_names(param_struct) 
 tagsize=size(tags)
 
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; create the command string
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 cmd_string = '~/SExtractor/mysextractor2013/source/sex '
 cmd_string = cmd_string+filename
 cmd_string = cmd_string + $
   ' -c ~/SExtractor/mysextractor2013/config/rotse.sex '

 n = n_elements(tags)
 ;; loop through the input tags and add them to command string
 for i = 0,n-1 do begin
  	tag = tags(i)
	val=string(param_struct.(i))
	cmd_string = cmd_string+'-'+tag+' '+val+' '
 endfor
   
 ;; send the command string to the shell
 spawn, cmd_string

 return

 end


