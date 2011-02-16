PRO read_corr, dir, struct, start=start, nframes=nframes, noadd=noadd, front=front, nchar=nchar, struct_type=struct_type, all=all, verbose=verbose

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: read_corr, dir, struct, start=start, nframes=nframes, noadd=noadd, front=front, nchar=nchar, struct_type=struct_type, all=all, verbose=verbose'
      return
  ENDIF 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time = systime(1)
  IF n_elements(verbose) EQ 0 THEN verbose=2
  IF n_elements(struct_type) EQ 0 THEN struct_type = 'sdss'
  IF NOT keyword_set(nomodrow) THEN nomodrow=0
  IF n_elements(noadd) EQ 0 THEN noadd=0
  IF n_elements(front) EQ 0 THEN front = 'adatc'
  IF n_elements(nchar) EQ 0 THEN nchar = 25

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Get filenames
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF keyword_set(all) THEN BEGIN
      fetch_file_list, dir, files, fnums, fieldmin=fieldmin, $
                       fieldmax=fieldmax, run=run, camcol=camcol, $
                       rerun=rerun, front=front, nchar=nchar
      startf = fieldmin
      nfields = fieldmax - startf + 1
  ENDIF ELSE BEGIN 
      IF n_elements(nframes) EQ 0 THEN nfields=1 ELSE nfields=nframes
      IF n_elements(start) EQ 0 THEN startf=0 ELSE startf=start
      fetch_file_list, dir, files, fnums, start=startf, nframes=nfields, $
                       fieldmin=fieldmin, fieldmax=fieldmax, $
                       run=run, camcol=camcol, rerun=rerun, $
                       front=front, nchar=nchar
  ENDELSE 

  IF files[0] EQ '' THEN return

  nfields = n_elements(files)
  runstr=strtrim(string(run), 2)
  colstr=strtrim(string(camcol),2)
  rerstr=strtrim(string(rerun),2)
  stop = ntostr(max(fnums))

  IF verbose GT 0 THEN BEGIN 
      print,' '
      print,'*****************************************'
      print,'* Run: ',runstr,'  Camcol: ',colstr,'  Rerun: ',rerstr
      print,'* Reading fields ',strtrim(string(startf),2),'-',stop
      print,'*****************************************'
      print,''
  ENDIF

  numlist=lonarr(nfields)
  ptrlist=ptrarr(nfields)       ;the array of pointers
  ntotal=0L
  
  FOR i=0L, nfields-1 DO BEGIN 
      
      openr, unit, files[i], /get_lun, /block, ERROR = error
      IF ERROR NE 0 THEN BEGIN 
          print,!ERR_STRING
          free_lun,unit
          return
      ENDIF 
      
      IF i EQ 0 THEN BEGIN
          lnew = mrdfits3(unit, 1, 0, /silent)
      ENDIF ELSE BEGIN
          lnew = mrdfits3(unit, 1, 0, /silent, /deja_vu)
      ENDELSE 
      free_lun, unit

      field = fnums[i]
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Due to overlap of consecutive images, throw out 
      ;; first 64 and last 64 rows unless first frame (i eq fieldmin)
      ;; or last frame (i eq fieldmax)
      ;; Can use nomodrow if reading files where this has already
      ;; been done!
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF NOT nomodrow THEN BEGIN 
          g=where(lnew.objc_rowc GT (field GT fieldmin)*64 AND $
                  lnew.objc_rowc LT (1489 - (field LT fieldmax)*64) , n)
      ENDIF ELSE BEGIN
          n=n_elements(lnew)
          g=lindgen( n )
      ENDELSE 

      ln = lnew[g]

      IF (NOT noadd) THEN lnew.objc_rowc = lnew.objc_rowc + $
                                             (1361L*(field-fieldmin))

      IF verbose GT 1 AND (  ((i MOD 20)  EQ 0 ) OR $
                            (i EQ nfields-1) OR (i EQ 0) ) THEN BEGIN 
          f=strtrim(string(field),2)
          o=strtrim(string(n),2)
          print,'Field: ',f,'/',stop,'  Objects: ',o
      ENDIF 

      numlist[i]=n
      ntotal=ntotal+n
      ptrlist[i]=ptr_new(ln)
                                ;put the pointer in the array of pointers
  ENDFOR 
   
  IF verbose GT 0 THEN BEGIN 
      o=strtrim(string(ntotal), 2)
      print,'Total Number of Objects: ',o
  ENDIF 
  struct=replicate(lnew[0],ntotal)
                                ;the output structure of correct size

;now make big list from the array of pointers ptrlist

  beg=0
  FOR i=0,nfields-1 DO BEGIN 
      struct(beg:beg+numlist(i)-1)=*ptrlist(i)
      ptr_free,ptrlist(i)       ;kill the heap variables associated 
                                ;with the pointer
      beg=beg+numlist(i)
  ENDFOR 

  ptime,systime(1)-time

return
END 
