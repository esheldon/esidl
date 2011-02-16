pro read_tsobjmin, file, field, fieldmin, fieldmax, indices, subl, $
                   struct, hdr, $
                   addrow=addrow,$
                   struct_type=struct_type, noadd=noadd, verbose=verbose, $
                   nomodrow=nomodrow, front=front, nchar=nchar

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
; NAME:  
;    read_tsobjmin
;
; PURPOSE:
;   minimal version of read_tsobj for reading one field at a time.  Should
;   be used with initcol.pro	
;
;    Notes: objc_rowc is made continuous from the first field rather than 
;      artificially restarting at the beginning of each field.  This can be 
;      easily removed by setting /noadd.  rowc[i] are not altered.  
; 
; CALLING SEQUENCE: 
;     read_tsobj, dir, struct, start=start,nframes=nframes,
;	taglist=taglist, phototags = phototags,addrow=addrow,
;	struct_type=struct_type,noadd=noadd, verbose=verbose,
;       front=front, nchar=nchar
;
; 
; INPUTS:  
;     file: file name (e.g. tsObj-*)
;     field: field number
;     fieldmin: first field
;     fieldmax: last field
;     indices: indices in of tags in subl to use
;     subl: structure containing desired tags
;     
;
; OPTIONAL INPUTS:
;
;          addrow:   Optional parameter which tells how many rows the user
;                    wants to add to objc_rowc (That tag must be in taglist)
;                    This may be necessary if there is more that one file
;                    for each column.
;	   struct_type: The user may specify a _different_ name for the
;			structure to be read in, if the file for example
;			does not contain default photo information
;          noadd=noadd: set this keyword if you don't want to add up 
;                       objc_rowc continuously.
;          verbose: Verbosity level. Default is verbose > 1, full.
;                       verbose = 0    print only error messages
;                       verbose = 1    print minimun of messages
;                       verbose > 1    print all messages.  Includes the 
;                                      field by field updates.
;       front:  The front string in the filename.  Default is "tsObj"
;       nchar:  The number of characters in the name.  
;               Default is 25 for tsObj files.
;
; OUTPUTS: struct:   The structure with desired tags
;
; OPTIONAL OUTPUTS: hdr: the header
;
; Example script:  Lets say you only want to process one field at a time.  That
;    is what this program is all about.
;
; run=752 & camcol=3 & rerun=1
; fetch_dir, run, camcol, rerun, dir
; 
; initcol, dir, files, fnums, indices, subl
; fieldmin = min(fnums) & fieldmax = max(fnums)
;
; start = fieldmin
; nframes = n_elements(fnums)
;
; for i=0, nframes-1 do begin
;
;    field = start+i
;    file  = files[field]
;    
;    read_tsobjmin, file, field, fieldmin, fieldmax, indices, subl, struct
;
;    (do stuff with struct....)
; endfor
; 
; Author:  Erin Scott Sheldon
; Date: 2/3/00
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Help message
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  IF n_params() EQ 0 THEN BEGIN 
      print,'-syntax read_tsobj, dir, struct, '
      print,'    start=start,nframes=nframes, all=all, '
      print,'    taglist=taglist, phototags = phototags,'
      print,'    addrow=addrow, struct_type=struct_type, noadd=noadd,'
      print,'    verbose=verbose, front=front, nchar=nchar'
      print,''
      print,'Use doc_library,"read_tsObj"  for more help'
      return
  ENDIF 

  time = systime(1)

  IF n_elements(verbose) EQ 0 THEN verbose=2
  IF NOT keyword_set(nomodrow) THEN nomodrow = 0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  read in this field
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF (NOT keyword_set(noadd) ) THEN noadd=0 ELSE noadd=1
  IF n_elements(struct_type) EQ 0 THEN struct_type = 'sdss'
	
  openr, unit, file, /get_lun, /block, ERROR = error
  IF ERROR NE 0 THEN BEGIN 
      print,!ERR_STRING
      free_lun,unit
      return
  ENDIF 

                                ;Must be setup previously!
  lnew = mrdfits3(unit,1,0,hdr,structyp=struct_type,$
                  /silent,/deja_vu)
                                ;it can use the structure that it already 
                                ;made the first time around
  free_lun,unit

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Due to overlap of consecutive images, throw out 
  ;; first 64 and last 64 rows unless first frame (i eq fieldmin)
  ;; or last frame (i eq fieldmax)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  IF NOT nomodrow THEN BEGIN 
      g=where(lnew.objc_rowc GT (field GT fieldmin)*64 AND $
              lnew.objc_rowc LT (1489 - (field LT fieldmax)*64),n)
  ENDIF ELSE BEGIN
      n=n_elements(lnew)
      g=lindgen( n )
  ENDELSE 
  IF (NOT noadd) THEN lnew.objc_rowc=lnew.objc_rowc+(1361L*(field-fieldmin))

  f=strtrim(string(field),2)
  o=strtrim(string(n),2)
  print,'Field: ',f,'  Objects: ',o

  struct = 0
  struct=replicate(subl,n)
  FOR k=0, n_elements(indices)-1 DO BEGIN 
                ;;;;must add 4 due to field, camcol, rerun, and run;;;;;
      struct.(k+4) = lnew(g).(indices[k])
  ENDFOR 
  struct.run = subl.run
  struct.camcol = subl.camcol
  struct.field = field
  struct.rerun = subl.rerun


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; May need to add additional rows since there are multiple 
; files for each column
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(addrow) THEN BEGIN 
      w = where(taglist eq "OBJC_ROWC")
      IF (w[0] NE -1) THEN BEGIN
          IF verbose GT 0 THEN BEGIN 
              print,'========================================='
              print,'Adding ',addrow,' rows to objc_rowc'
              struct.objc_rowc = struct.objc_rowc + addrow
          ENDIF 
      ENDIF 
  ENDIF 

  return
END 
	






