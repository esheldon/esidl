PRO set_bbox, filename, bbox_string
  
  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: set_bbox, filename, bbox_string'
      print
      print,'e.g. set_bbox, psfile, "%%BoundingBox: 35 15 770 515"'
      return
  ENDIF 

  tempfile = tmpfile(prefix='/tmp/tmpfile_')

  openr, lun,    filename, /get_lun
  openw, tmplun, tempfile,  /get_lun

  lin=''
  replinz=strarr(100)
  linno=lonarr(100)
  ctr=0 & linctr=0L

  search_string = '%%BoundingBox'

  while not EOF(lun) do begin
      readf,lun,lin
      mm = strpos(lin, search_string)
      if mm[0] NE -1 then begin
          print,lin+' -> '+bbox_string
          printf, tmplun, bbox_string
          print,'---'

      ENDIF ELSE BEGIN 
          printf, tmplun, lin
      ENDELSE 
  endwhile
  free_lun,lun
  free_lun, tmplun

  spawn,'cp -f '+tempfile+' '+filename+'; rm -f '+tempfile

END 
