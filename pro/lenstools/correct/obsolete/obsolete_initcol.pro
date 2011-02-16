PRO initcol, dir, files, fnums, indices, subl, lensstr, $
             start=start, nframes=nframes, $
             fieldmin=fieldmin, fieldmax=fieldmax

  IF n_params() LT 1 THEN BEGIN
      print,'-Syntax: initcol, dir, files, fnums, indices, subl, lensstr, '
      print, '         start=start, nframes=nframes, '
      print,'         fieldmin=fieldmin, fieldmax=fieldmax'
      return
  ENDIF 

  fetch_file_list, dir, files, fnums, run=run, camcol=camcol, rerun=rerun, $
                   start=start, nframes=nframes, $
                   fieldmin=fieldmin, fieldmax=fieldmax

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure out which taglist to use
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  make_lenstags,taglist,/default
  ;; add these just for finding E_D_BIT
  taglist = [ taglist, ['EXP_L', 'DEV_L'] ]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Get all the phototags and prep deja_vu
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  openr, unit, files(0), /get_lun, /block, ERROR = error
  IF ERROR NE 0 THEN BEGIN 
      print,!ERR_STRING
      free_lun,unit
      return
  ENDIF 
  
  lnew = mrdfits3(unit,1,0,hhdr,structyp=struct_type,/silent)
  free_lun,unit
  phototags = tag_names(lnew)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;figure out which tags are good and use them to make the struct
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  read_tsobj_make_tags,taglist,phototags,goodtags,indices
  IF goodtags[0] EQ '' THEN return
  
  read_tsobj_make_struct, goodtags, indices, lnew, subl

  subl.run = run
  subl.camcol = camcol
  subl.rerun = rerun

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Stuff to keep for lenstags
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  make_lenstags, taglist

  read_tsobj_make_tags,taglist,phototags,n_goodtags,n_indices
  IF goodtags[0] EQ '' THEN return
  
  read_tsobj_make_struct, n_goodtags, n_indices, lnew, n_subl

  n_subl.run = run
  n_subl.camcol = camcol
  n_subl.rerun = rerun

  lensstr = create_struct(n_subl, $
                          'STARFLAG', intarr(5), $
                          'E_D_BIT', 0, $
                          'CLASSIFICATION',intarr(5), $
                          'IXX',fltarr(5), $
                          'IYY',fltarr(5), $
                          'IXY',fltarr(5), $
                          'RHO4',fltarr(5), $
                          'WHYFLAG',intarr(5), $
                          'E1',fltarr(5), $
                          'E2',fltarr(5), $
                          'MOMERR',fltarr(5), $
                          'R',fltarr(5), $
                          'PHOTOZ', 0. )

return
END 
