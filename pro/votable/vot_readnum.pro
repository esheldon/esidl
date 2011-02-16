;+
; NAME:
;  VOT_READNUM()
;
;
; PURPOSE:
;  Read a VO table containing only numbers and arrays of numbers into an IDL
;  structure.   Because it is only numbers, this can be done very efficiently.
;  Currently only supports a single table in the file, and assumes the field
;  definitions are each contained on a single line, not including the
;  descripton.  See restrictions below.
;
;
; CATEGORY:
;  Virtual Observatory XML table reader.
;
; CALLING SEQUENCE:
;    struct = vot_readnum(file)
;
; INPUTS:
;  file: The vo table file.
;
; OUTPUTS:
;  An idl structure is returned.
;
; OPTIONAL OUTPUTS:
;  status=:  status !=0 means execution failed
;
; KEYWORDS:
;  /silent:  Don't print informative messages.
;
; COMMON BLOCKS:
;  vot_readnum_COMMON, fstring, fstring_len, tabstring, tabstring_len
;
; RESTRICTIONS:
;  1) You must have compiled the DLM remove_xml_tags and it must be in
;     your IDL_DLM_PATH
;  2) Only a single table per file is currently valid.
;  3) Field definitions must be on a single line. e.g.
;    <FIELD arraysize="5" datatype="float" name="petroCounts"/>
;  or
;    <FIELD name="vmag" datatype="float">
;
; EXAMPLE:
;  file = 'tsObj-000756-6-44-0367.vot'
;  struct = vot_readnum(file)
;
;
; MODIFICATION HISTORY:
;  Created: 4-March-2005  Erin Sheldon, UChicago
;
;-


PRO vot_readnum_com

  COMMON vot_readnum_COMMON, fstring, fstring_len, tabstring, tabstring_len

  fstring='<field'
  fstring_len = 6
  tabstring='>'
  tabstring_len = 1

END 


PRO vot_readnum_parse_field_def, field_def, tagname, tagtype

  COMMON vot_readnum_COMMON, fstring, fstring_len, tabstring, tabstring_len

  t=strlowcase(field_def)

  ;; Extract from between the <field ... >
  p1 = strpos(t, fstring)
  p2 = strpos(t, tabstring,/reverse_search)

  sstart = p1 + fstring_len
  slen = p2-sstart

  t = strmid(temporary(t), sstart, slen)

  ;; Split on white space
  t=strsplit(temporary(t), /extract)

  ;; Get the tag name
  w=where(strmatch(t,"name*"),nw)
  tagname = strsplit(t[w[0]], '"', /extract)
  tagname = tagname[1]

  ;; Get the array size if there and the 
  ;; tag type

  w=where(strmatch(t, "arraysize*"), nw)

  IF nw EQ 0 THEN BEGIN 

      ;; Scalar
      w=where(strmatch(t, "datatype*"), nw)

      tp = t[w[0]]
      tp = strsplit(tp, '"', /extract)
      tp = tp[1]

      CASE tp OF
          'short':  tagtype='0'
          'int':    tagtype='0L'
          'long':   tagtype='0LL'
          'float':  tagtype='0.0'
          'double': tagtype='0d'
          ELSE: message,'Unsupported data type: '+tp
      ENDCASE 

  ENDIF ELSE BEGIN 

      ;; Array size
      arrsz = t[w[0]]
      arrsz = strsplit(arrsz, '"', /extract)
      arrsz = arrsz[1]

      arrsz = repstr(arrsz, "x", ",")

      w=where(strmatch(t, "datatype*"), nw)

      tp = t[w[0]]
      tp = strsplit(tp, '"', /extract)
      tp = tp[1]

      CASE tp OF
          'short':  tagtype='intarr('+arrsz+')'
          'int':    tagtype='lonarr('+arrsz+')'
          'long':   tagtype='lon64arr('+arrsz+')'
          'float':  tagtype='fltarr('+arrsz+')'
          'double': tagtype='dblarr('+arrsz+')'
          ELSE: message,'Unsupported data type: '+tp
      ENDCASE 
  ENDELSE 
  
  
END 

FUNCTION vot_readnum_getstruct, lines, silent=silent

  ;; Get the data types from the field definitions
  ;; and create a structure

  nl=n_elements(lines)
  tag_names = strarr(nl)
  tag_types = strarr(nl)
  FOR i=0L, nl-1 DO BEGIN 
      vot_readnum_parse_field_def, lines[i], tagname, tagtype
      tag_names[i] = tagname
      tag_types[i] = tagtype
  ENDFOR 

  struct = mrd_struct(tag_names, tag_types, 1)
  return,struct
END 

FUNCTION vot_readnum_openfile, file, error=error

  openr, lun, file, /get_lun, error=error

  IF error NE 0 THEN BEGIN 
      print,'Error opening file '+file+': '+!error_state.sys_msg
      return,-1
  ENDIF ELSE BEGIN 
      return,lun
  ENDELSE 

END 

FUNCTION vot_readnum, file, silent=silent, status=status

  status=1

  COMMON vot_readnum_COMMON, fstring, fstring_len, tabstring, tabstring_len
  IF n_elements(fstring) EQ 0 THEN vot_readnum_com

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: struct = vot_readnum(file, /silent, status=)'
      return,-1
  ENDIF 

  ;; Read the file into an array of strings
  IF NOT keyword_set(silent) THEN print,'Reading file: ',file

  nlines = numlines(file)
  lun = vot_readnum_openfile(file, error=error)
  IF error NE 0 THEN return,-1

  lines = strarr(nlines)
  readf, lun, lines
  free_lun, lun

  ;; Get the field definitions
  nlines = n_elements(lines)
  i=0L
  ibeg = -1L
  WHILE 1 DO BEGIN 

      IF strmatch(lines[i], "<FIELD*", /FOLD_CASE) THEN BEGIN 
          add_arrval, i, fieldi
      ENDIF ELSE BEGIN 
          IF strmatch(lines[i], "<TABLEDATA*",/FOLD_CASE) THEN BEGIN 
              tabi = i
              BREAK 
          ENDIF 

      ENDELSE 
      
      i=i+1
  ENDWHILE 

  ;; Use the field descriptions to create an IDL struct
  st = vot_readnum_getstruct(lines[fieldi])

  ;; Grab the substring. This is kind of a bottleneck
  lines = lines[tabi:nlines-1]

  ;; strjoin is very efficient
  lcat = strjoin(lines)

  ;; This is by far the fastest way to find the endpoints
  ;; Only works when there is a single table
  s1 = '<TABLEDATA>'
  s1len = strlen(s1)
  t1 = strpos(lcat, s1)


  s2 = '</TABLEDATA>'
  s2len = strlen(s2)
  t2 = strpos(lcat, s2, /reverse_search)

  ;; extract the substring
  sstart = t1 + s1len
  slen = t2-sstart
  lcat = strmid(temporary(lcat), sstart, slen)


  ;; Now we have a bunch of rows with <TR> .... </TR>
  ;; and within them we have a bunch of datas <TD> .... </TD>
  ;; We want to remove these tag definitions from the string entirely
  ;; in order to use the fast IDL procedure reads to read from the string
  ;; directly into a structure
  
  ;; remove_xml_tags is a DLM written in C that removes all tags <..>
  ;; and counts the <TR>s, or number of rows, along the way
  ;; Note, the temporary saves memory becuase the operation is done
  ;; in-place in the original string

  data = remove_xml_tags(temporary(lcat), nrows=nrows)

  IF NOT keyword_set(silent) THEN print,'nrows = ',nrows

  struct = replicate(st, nrows)
  reads, data, struct

  ;; free the memory
  data=0

  status = 0
  return,struct

END 
