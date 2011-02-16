PRO newcorr, run, camcol, rerun, $
             newfront=newfront, oldfront=oldfront

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: addsdsspar, run, camcol, rerun, newfront=newfront, oldfront=oldfront'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; check params
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(oldfront) EQ 0 THEN oldfront = 'adatc'
  IF n_elements(newfront) EQ 0 THEN newfront = 'adatcadd' 

  print
  print,'Old front: ',oldfront
  print,'New front: ',newfront
  print

  fetch_dir, run, camcol, rerun, dir, corrdir=corrdir

  oldnchar = 20 + strlen(oldfront)
  newnchar = 20 + strlen(newfront)

  fetch_file_list, dir, tsfiles, tsfnums
  fetch_file_list, corrdir, corrfiles, fnums, front = oldfront, nchar=oldnchar

  rename_tsobj, tsfiles, corrdir, newfront, newfiles

  initcol, dir, ff, fnn, iinndd, sl, lensstr

  taglist = tag_names(sl)

  nf = n_elements(fnums)

  time=systime(1)

  FOR i=0, nf-1 DO BEGIN 
      

      field = fnums[i]
      outfile = newfiles[i]

      struct=0
      read_tsobj, dir, struct, start=field, taglist=taglist, tsobjstr=tsobjstr

      corstruct = 0
      corrstruct = mrdfits(corrfiles[i], 1, hdr, /silent)
      FXHCLEAN,hdr

      ncorr = n_elements(corrstruct)
      
      photo_match, struct.run, struct.rerun, struct.camcol, struct.field, $
        struct.id, corrstruct.run, corrstruct.rerun, corrstruct.camcol, $
        corrstruct.field, corrstruct.id, $
        match, corrmatch

      ncorrmatch = n_elements(corrmatch)
      IF ncorrmatch NE ncorr THEN BEGIN
          print,'WHat!'
          return
      ENDIF 

      newstruct = 0
      newstruct = replicate(lensstr, ncorr)

      copy_struct, struct[match], newstruct
      copy_struct, corrstruct[corrmatch], newstruct

      ;; This part may be removed after this run through?
      newstruct.e_d_bit = 1
      wdev = where(struct[match].dev_l[2] GT struct[match].exp_l[2],ndev)
      IF ndev NE 0 THEN newstruct[wdev].e_d_bit = 2

      diff_list = compare_struct(corrstruct[corrmatch],newstruct)
      IF diff_list.ndiff NE 0 THEN BEGIN
          help,diff_list,/str
          return
      ENDIF 

;      help,newstruct[wdev[0]],/str

      print,outfile
      mwrfits, newstruct, outfile, hdr, /create
;      return

      
  ENDFOR 
  print,'----------------'
  ptime,systime(1)-time

  return
END 
