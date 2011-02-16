;+
; NAME:
;  WRITE_IDLSTRUCT
;
;
; PURPOSE: 
;   Write a structure array to a standardized file type ".st".  The file is
;   self-describing such that it can be read using READ_IDLSTRUCT() without
;   prior knowledge of the file structure.  The data is either written as
;   unformatted binary in the native or IEEE formats or as ascii when the
;   /ascii keyword is sent.  These procedures provide a natural way to write
;   and read IDL data in ascii or binary form.  Unlike fits files, all IDL data
;   types other than complex, pointers, and sub-structures are supported for
;   both formats, including array fields for strings.  Data can be appended
;   naturally using the /append keyword.
;
;
; CATEGORY:
;  File I/O
;
; CALLING SEQUENCE:
;  write_idlstruct, struct, filename/unit, /ascii, /ieee, /append, 
;                   hdrStruct=hdrStruct, 
;                   info_struct=info_struct
;
;
; INPUTS:
;
;  struct: An IDL structure array. All basic IDL data types other than
;          complex, pointers, and sub-structures are supported for both binary
;          and ascii formats.
;
;  filename/unit: Name of the file to be written or appended, or the logical
;                 unit of an open file. Array fields in the structure is
;                 supported except for strings.
;
; OPTIONAL INPUTS:
;  hdrStruct=: The user can input keyword-value to be written into the header
;              through this structure.  For example:
;                   hdr={date:"12-May-1974", val1:[2.5, 66.7], val2:33}
;                   write_idlstruct, struct, file, hdrStruct=hdr
;
;  info_struct=: This structure contains info on the makeup of the structure,
;               and output formatting for ascii.  It can be returned on the
;               first call and sent again when appending to save time. Even
;               with this time saver, writing line-by-line is very slow.
;
;
; KEYWORD PARAMETERS:
;  /ascii: write an ascii data file rather than unformatted binary.
;  /ieee: Use IEEE XDR byte ordering rather than the native for
;           LITTLE_ENDIAN machines (they are the same for BIG_ENDIAN) Using
;           native format can increase the speed of reads and writes by a
;           large amount, but then BIG_ENDIAN users must swap the byteorder to
;           use the data.  Also, read_idlstruct will convert to IEEE if needed,
;           so this is not necessary except for speed on the BIG_ENDIAN
;           machines. 
;  /append: append to the file and update the header.  WARNING: currently,
;           appending will not work with ascii files that contain strings.
;           A fixed field width must be used, and this procedure does not
;           currently check that when appending.
;
; OUTPUTS:
;  The structure is written to a file
;
;
; OPTIONAL OUTPUTS:
;  info_struct: see above
;
; PROCEDURES CALLE:
;  DATATYPE
;  (HOST_TO_IEEE)
;  (IEEE_TO_HOST)
;
; SIDE EFFECTS:
;  For string fields in the structure, the strings are padded to the largest
;  string in the array for that field.
;
; COMMON BLOCKS:
;  write_idlstruct_block, ROW_STRING_FORMAT
;
; EXAMPLE:
;  struct = {a:findgen(10), b:dindgen(10,10), c:23353233}
;  struct = replicate(struct, 1000)
;  write_idlstruct, struct, 'test.st'
;  write_idlstruct, struct, 'test2.st', /ascii
;
;
; MODIFICATION HISTORY:
;  Created  02-Jul-2004 Erin Sheldon, UofChicago
;
;-


FUNCTION write_postgres_input_openfile, filename, append=append, error=error, $
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

FUNCTION write_postgres_input_format, struct, noarrays=noarrays, error=error

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
                  ;; oformat = '%d'
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

PRO write_postgres_input, struct, filename, $
                          format=format, append_in=append_in, $
                          noarrays=noarrays, $
                          error=error

  error = -20000

  IF n_params() LT 2 THEN BEGIN 
    print,'-Syntax: write_postgres_input, struct, filename/unit, $'
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
  lun = write_postgres_input_openfile(filename, append=append, error=error, $
                                 tname=tname)
  IF error NE 0 THEN return

  ;; If the info struct not sent, get the info we need
  ;; Otherwise, just update the number of rows

  nrows = n_elements(struct)

  IF n_elements(format) EQ 0 THEN BEGIN 
      format = $
        write_postgres_input_format(struct,noarrays=noarrays,error=error)
      IF error NE 0 THEN return
  ENDIF 

  ;; Output the structure
  FOR row=0ULL, nrows-1 DO BEGIN 
      printf,lun,struct[row],format=format
  ENDFOR 

  IF (tname EQ 'STRING') THEN free_lun,lun

END 
