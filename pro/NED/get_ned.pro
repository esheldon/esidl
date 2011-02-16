PRO get_ned, ned, typestring, typestruct, file=file, plot=plot

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: 
;    GET_NED
;       
; PURPOSE: 
;    read in the ned catalog from the specified FITS file.  
;    NOTE: if you need to read it in from the ascii file, use
;          read_ned_data.  Then output it as a fits file and 
;          use this reader.  Much faster that ascii reading.
;	
;
; CALLING SEQUENCE:
;      get_ned, ned, other, otherstring, file=file, plot=plot
;                 
;
; OPTIONAL INPUTS:
;      typestring: Name of special type that you want put into structure 
;                  typestr.
;      file: the ned fits file
;  
; KEYWORD PARAMETERS:
;      /plot: plot ra's and decs.
;       
; OPTIONAL OUTPUTS: ned: the ned catalog
;          typestruct: catalog of special type objects.
;
; REVISION HISTORY:
;   Author:  Erin Scott Sheldon  5/??/99	
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() EQ 0 THEN BEGIN 
     print,'-Syntax: get_ned, [, ned, typestring, typestruct, file=file, plot=plot]'
     print,''
     print,'Use doc_library,"get_ned"  for more help.'  
     return
  ENDIF 

  IF n_elements(typestring) NE 0 THEN gettype = 1 ELSE gettype=0
  IF NOT keyword_set(file) THEN BEGIN
    file='/sdss3/data2/Ned/sdss_equatorial_with_z.fit'
  ENDIF 

  ned = mrdfits(file,1,hdr,/silent)
  n=n_elements(ned)
  FOR i=0,n-1 DO BEGIN
    ned[i].type = strtrim(ned[i].type)
    ned[i].name = strtrim(ned[i].name)
  ENDFOR
  help,ned
  print,gettype
  s = ned[ rem_dup(ned.type) ].type

  IF gettype THEN BEGIN
    w=where(ned.type EQ typestring, nw)
    IF nw EQ 0 THEN BEGIN
        print,'Type: '+typestring+' Not Found'
        print,'Try one of these: ',s
        gettype = 0
    ENDIF ELSE BEGIN 
        typestruct = ned[w]
        help,typestruct
    ENDELSE 
  ENDIF

  IF keyword_set(plot) THEN BEGIN
    plot,ned.ra,ned.dec,psym=3,xtitle='ra',ytitle='dec',title='Ned'
    IF gettype THEN oplot,typestruct.ra,typestruct.dec,psym=7
  ENDIF
  
  return
end



