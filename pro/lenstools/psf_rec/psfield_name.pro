pro psfield_name, run, camcol, field, fname

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: 
;    PSFIELD_NAME
;
; PURPOSE: 
;    
;
; CALLING SEQUENCE:
;    
;
; INPUTS:  
;    run: sdss run number
;    camcol: camera column
;    field: field of the object
;
; OUTPUTS: 
;    psfieldname
;
; REVISION HISTORY:
;    Erin Scott Sheldon  02/23/00
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: psfield_name, run, camcol, field [, fname]'
      return
  ENDIF 

  rs = run2string(run)
  cs = strtrim(string(camcol),2)
  fs = field2string(field)
  fname='psField-'+rs+'-'+cs+'-'+fs+'.fit '
  fname = strtrim(fname)

return
end
