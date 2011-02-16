pro  read_rosat, file, str

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: read_rosat
;       
; PURPOSE: read contents of input file into a structure.
;	
;
; CALLING SEQUENCE: read_rosat, file, str
;
; INPUTS: file:   file name string.
;       
; OUTPUTS: str:  ouput structure containing rosat info.
; 
; PROCEDURE: 
;	Use predefined structure (gotten from file) and read file into it.
;	
;
; REVISION HISTORY:
;	Erin Scott Sheldon   Umich 5/25/99
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if N_params() eq 0 then begin
     print,'-Syntax: read_rosat, file, str'
     print,''
     print,'Use doc_library,"read_rosat"  for more help.'  
     return
  endif



get_lun, flun
openr, flun, file

string=''

s = create_struct('srcnam', '', $
                  'rosat_ra', 0.0, $
                  'rosat_dec',0.0, $
                  'pos_err',0, $    ;an integer, includes 6" systematic error
                  'npedm_flags', '', $
                  'riv_flags', '', $
                  'src_cps', 0.0, $
                  'src_cps_err',0.0, $
                  'bgrcpsa',0.0, $
                  'exp', 0L, $
                  'hr1',0.0, $
                  'hr1_err',0.0, $
                  'hr2', 0.0, $
                  'hr2_err', 0.0, $
                  'ext', 0.0, $
                  'ext1',0.0, $
                  'src1', 0.0, $
                  'extr', 0.0, $
                  'PriFlg', 0L, $
                  'E', '', $
                  'vigf', 0.0, $
                  'orgdat', '', $
                  'moddat', '', $
                  'field_id', 0L, $
                  'src_num', 0L)
                  
nlines = 18811
;nlines=33
ncolumns = '185'
format = '(A' + ncolumns + ')'

str = replicate(s, nlines)

FOR i=0L, nlines-1 DO BEGIN
  readf, flun, string, format=format

  str[i].srcnam = strmid(string, 0, 21)
  str[i].rosat_ra = double( strmid(string,23,9) )
  str[i].rosat_dec = double( strmid(string,33,9) )
  str[i].pos_err = long( strmid(string,44, 2) )

  str[i].npedm_flags = strmid(string,47,5)

  str[i].riv_flags = strmid(string,53,3)

  str[i].src_cps = float(strmid(string,59, 9))
  str[i].src_cps_err = float(strmid(string,69, 9))
  str[i].bgrcpsa = float(strmid(string,79,9))
  str[i].exp = long(strmid(string, 89,6))
  str[i].hr1 = float(strmid(string,96,5))
  str[i].hr1_err = float(strmid(string,102,4))
  str[i].hr2 = float(strmid(string,107,5))
  str[i].hr2_err = float(strmid(string,113,4))
  str[i].ext = float(strmid(string,118,5))
  str[i].ext1 = float(strmid(string,124,4))
  str[i].src1 = float(strmid(string, 129, 4))
  str[i].extr = float(strmid(string,134,5))
  str[i].PriFlg = long(strmid( string,140,6))
  str[i].E = strmid(string, 146,1)
  str[i].vigf = float( strmid(string,148,4) )
  str[i].orgdat = strmid(string,153,6)
  str[i].moddat = strmid(string,160,6)
  str[i].field_id = long( strmid(string, 172,8))
  str[i].src_num = long( strmid(string, 181, 4))

  x = fix(i/1000) - i/1000.0 
  IF ( x EQ 0 ) THEN BEGIN
    print,'Have read ',strtrim(string(i), 2)
  ENDIF 

ENDFOR 


return
end








