;+
; NAME:
;  POSTGRES_WRITE_INPUT
;
;
; PURPOSE: 
;   Write a structure array to a file suitable for import into PostgreSQL.
;   Columns are tab-separated. 
;
; CATEGORY:
;  File I/O
;
; CALLING SEQUENCE:
;  postgres_write_input, struct, filename, format=, append=, $
;                        /noarrays, error=
;
; INPUTS:
;
;  struct: An IDL structure array. All basic IDL data types other than
;          complex, pointers, and sub-structures are supported.
;
;  filename/unit: Name of the file to be written or appended, or the logical
;                 unit of an open file. 
;
; OPTIONAL INPUTS/OUTPUTS:
;  format=: the format string for printing to the file.  By default, this is
;           generated within the program. It can be returned and send again
;           during appends to speed things  up.
;  error=: Error code.  If error=0, execution was successful.
;
; KEYWORD PARAMETERS:
;  /append: append to the file and update the header. 
;  /no_arrays: By default, arrays are written so that they will be stored in
;              postgres as an array.  E.g. [1,2,3] -> {1,2,3} in the file. If
;              this keyword is set, the are written as normal columns. 
;
; EXAMPLE:
;  struct = {a:findgen(10), b:dindgen(10,10), c:23353233}
;  struct = replicate(struct, 1000)
;  postgres_write_input, struct, 'test.st'
;
; MODIFICATION HISTORY:
;  Created  01-May-2005 Erin Sheldon, UofChicago
;
;-


FUNCTION postgres_write_input_openfile, filename, append=append, error=error, $
                                   tname=tname

  tname = size(filename,/tname)
  IF (tname EQ 'INT' OR tname EQ 'UINT' OR $
      tname EQ 'LONG' OR tname EQ 'ULONG' OR $
      tname EQ 'LONG64' OR tname EQ 'ULONG64' OR $
      tname EQ 'BYTE') THEN BEGIN 

      lun = filename
      error = 0
      return, lun

  ENDIF ELSE IF tname EQ 'STRING' THEN BEGIN 

      IF keyword_set(append) AND NOT fexist(filename) THEN append=0
      openw, lun, filename, /get_lun, error=error, append=append

      IF error NE 0 THEN BEGIN 
          print,'Error opening file '+filename+': '+!error_state.sys_msg
          return,-1
      ENDIF ELSE BEGIN 
          return,lun
      ENDELSE 

  ENDIF ELSE BEGIN
 
      print,'Datatype of filename/unit '+tname+' is incorrect'
      error = -2000
      return,-1

  ENDELSE 

END 

FUNCTION postgres_write_input_format, struct, noarrays=noarrays, error=error

  tags = strlowcase( tag_names(struct) )
  ntags = n_elements(tags)

  tmp = struct[0]
  zero_struct, tmp

  tab = '"'+string(9b)+'"'

  format = '('
  FOR tag=0L, ntags-1 DO BEGIN 
      tstr = size(tmp.(tag), /struct)
      type = tstr.type_name
      desc = datatype(tmp.(tag),/desc)

      nel = tstr.n_elements

      IF NOT keyword_set(noarrays) AND nel GT 1 THEN format = format + '"{",'

      ;; For arrays
      FOR j=0L, nel-1 DO BEGIN 
          CASE type OF
              'INT': BEGIN 
                  format = format + 'I0'
              END 
              'LONG': BEGIN 
                  format = format + 'I0'
              END 
              'UINT': BEGIN 
                  format = format + 'I0'
              END 
              'ULONG': BEGIN 
                  format = format + 'I0'
              END 
              'LONG64': BEGIN 
                  format = format + 'I0'
              END 
              'ULONG64': BEGIN 
                  format = format + 'I0'
              END 

              'BYTE': BEGIN 
                  format = format + 'I0'
              END 

              'FLOAT': BEGIN 
                  format = format + 'g0'
              END 
              'DOUBLE': BEGIN 
                  format = format + 'e15.8'
              END 
              'STRING': BEGIN 
                  format = format + 'a'
              END 
              'COMPLEX': BEGIN 
                  print,'Complex type not supported'
                  error=-3000
                  return,-1
              END 
              ELSE: BEGIN 
                  print,'Type not supported: '+type
                  error = -4000
                  return,-1
              END 
          ENDCASE 
          ;; arrays
          IF nel GT 1 AND j NE nel-1 THEN BEGIN 
              IF NOT keyword_set(noarrays) THEN format=format+',",",' $
              ELSE format = format + ','+tab+','
          ENDIF 
      ENDFOR 

      IF NOT keyword_set(noarrays) AND nel GT 1 THEN format = format + ',"}"'
      ;; csv
;      IF tag NE ntags-1 THEN format = format + ',",",'

      ;; tab
      IF tag NE ntags-1 THEN format = format + ','+tab+','

  ENDFOR 

  format = format +')'

  return,format

END 

FUNCTION postgres_write_input_formats, struct, error=error

  tags = strlowcase( tag_names(struct) )
  ntags = n_elements(tags)

  formats = strarr(ntags)

  tmp = struct[0]
  zero_struct, tmp

  tab = '"'+string(9b)+'"'


  FOR tag=0L, ntags-1 DO BEGIN 
      tstr = size(tmp.(tag), /struct)
      type = tstr.type_name
      desc = datatype(tmp.(tag),/desc)

      nel = tstr.n_elements

      ;; No newline except on last tag
      IF tag NE ntags-1 THEN BEGIN 
          format = '($, '
      ENDIF ELSE BEGIN 
          format = '('
      ENDELSE 

      IF nel GT 1 THEN format = format + '"{",'

      ;; For arrays
      FOR j=0L, nel-1 DO BEGIN 
          CASE type OF
              'INT': BEGIN 
                  format = format + 'I0'
              END 
              'LONG': BEGIN 
                  format = format + 'I0'
              END 
              'UINT': BEGIN 
                  format = format + 'I0'
              END 
              'ULONG': BEGIN 
                  format = format + 'I0'
              END 
              'LONG64': BEGIN 
                  format = format + 'I0'
              END 
              'ULONG64': BEGIN 
                  format = format + 'I0'
              END 

              'BYTE': BEGIN 
                  format = format + 'I0'
              END 

              'FLOAT': BEGIN 
                  format = format + 'g0'
              END 
              'DOUBLE': BEGIN 
                  format = format + 'e15.8'
              END 
              'STRING': BEGIN 
                  format = format + 'a'
              END 
              'COMPLEX': BEGIN 
                  print,'Complex type not supported'
                  error=-3000
                  return,-1
              END 
              ELSE: BEGIN 
                  print,'Type not supported: '+type
                  error = -4000
                  return,-1
              END 
          ENDCASE 
          ;; arrays
          IF nel GT 1 AND j NE nel-1 THEN BEGIN 
              format=format+',",",' 
          ENDIF 
      ENDFOR 

      IF nel GT 1 THEN format = format + ',"}"'


      ;; Add a tab except for the last one
      IF tag NE ntags-1 THEN BEGIN 
          format = format + ','+tab
      ENDIF 

      format = format + ')'

      formats[tag] = format

  ENDFOR 

  return,formats

END 



PRO postgres_write_input, struct, filename, $
                          format=format, append_in=append_in, $
                          noarrays=noarrays, $
                          error=error

  error = -20000

  IF n_params() LT 2 THEN BEGIN 
    print,'-Syntax: postgres_write_input, struct, filename/unit, $'
    print,'                format=, /noarrays, error='
    return 
  ENDIF 

  ;; Copy this, since if file doesn't exist, then append will be set to 0
  IF n_elements(append_in) NE 0 THEN append=append_in
    
  ;; Check the structure
  IF datatype(struct, /tname) NE 'STRUCT' THEN BEGIN 
      print,'Data must be in structure form'
      error = -10000
      return
  ENDIF 

  ;; Open the file or just set lun if filename is file unit
  lun = postgres_write_input_openfile(filename, append=append, error=error, $
                                 tname=tname)
  IF error NE 0 THEN return

  ;; If the info struct not sent, get the info we need
  ;; Otherwise, just update the number of rows

  nrows = n_elements(struct)

  IF n_elements(formats) EQ 0 THEN BEGIN 
      formats = postgres_write_input_formats(struct,error=error)
      IF error NE 0 THEN return
  ENDIF 

;  colprint,formats

  ;; Output the structure
  nrows = n_elements(struct)
  ntags = n_elements(formats)
  FOR row=0LL, nrows-1 DO BEGIN 

      FOR tag=0L, ntags-1 DO BEGIN 
          printf, lun, struct[row].(tag), format=formats[tag]
      ENDFOR 

  ENDFOR 

  IF (tname EQ 'STRING') THEN free_lun,lun

  ;; This is faster, and works the same as long as no : are used
  ;; !!!! But doesn't work for long format strings ARGGG!!!!
;  printf,lun,struct,format=format



END 
