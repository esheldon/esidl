pro get_atlas_min, id, color_index, fname, im, s, row0, col0, maxsize=maxsize

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; NAME:
;       GET_ATLAS3
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if N_params() eq 0 then begin
      print,'-Syntax: get_atlas_min, id, color_index, fname, im, s, row0, col0, maxsize=maxsize'
      return
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Default size
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(maxsize) EQ 0 THEN size=[500,500] ELSE size=maxsize
  size_im = size
 
  defsysv, '!SOFILE', exists=exists1
  defsysv, '!ENTRY', exists=exists2
  IF (NOT exists1) OR (NOT exists2) THEN sdssidl_setup
  
  sofile = !SOFILE
  entry  = !ENTRY

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ids=strtrim(id,2)
  coli = strtrim(string(color_index),2)+' '
  
  im = intarr(size[0],size[1])
  origin=intarr(2)
 
  s=call_external(value=[0b,0b,0b,0b,0b,0b],sofile,entry,im,size_im,$
                  coli,fname,ids,origin)
  row0=origin[0]
  col0=origin[1]

  if (s eq 1) then begin
      return
  endif  

  im = temporary( im[ 0:size_im[0]-1, 0:size_im[1]-1 ] )

  return
end









