FUNCTION getkenturl, run, camcol, field, clusters=clusters

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    
;       
; PURPOSE:
;    
;
; CALLING SEQUENCE:
;    
;
; INPUTS: 
;    
;
; OPTIONAL INPUTS:
;    
;
; KEYWORD PARAMETERS:
;    
;       
; OUTPUTS: 
;    
;
; OPTIONAL OUTPUTS:
;    
;
; CALLED ROUTINES:
;    
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() EQ 0 THEN BEGIN 
     print,'-Syntax: url = getkenturl(run, camcol, field, clusters=clusters)'
     print,''
     print,'Use doc_library,""  for more help.'  
     return,''
  ENDIF 

  rstr = ntostr(long(run))
  cstr = ntostr(long(camcol))
  fstr = field2string(long(field))

  a = long(field)/10 & b= long(field) MOD 10
  IF b EQ 0 THEN BEGIN
      lowfield = (a-1)*10 + 1
      highfield = a*10
  ENDIF ELSE BEGIN
      lowfield = a*10 + 1
      highfield = (a+1)*10
  ENDELSE 
  
  frange = ntostr(lowfield)+'-'+ntostr(highfield)
  url = 'sdsslnx.fnal.gov:8015/'
  
  IF NOT keyword_set(clusters) THEN BEGIN 
      url = url + 'run-' + rstr + '/mosaic-' + frange 
      url = url + '/' + rstr + '/' + cstr
      url = url + '/h_' + fstr + '.htm'
  ENDIF ELSE BEGIN
      url = url + 'template/tsCluster.tml?root='
      url = url + '/run-'+rstr+'/mosaic-'+frange+'/'+rstr+'/'
      url = url + cstr+'&run='+rstr+'&camcol='+cstr
      url = url + '&field='+ntostr(long(field))
  ENDELSE 

  return, url
END 
