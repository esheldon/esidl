PRO addsdsspar, run, camcol, rerun, parlist, addstruct=addstruct, $
                newfront=newfront, oldfront=oldfront

  ;; adds parameter to corrected files

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: addsdsspar, run, camcol, rerun, parlist, addstruct=addstruct, newfront=newfront, oldfront=oldfront'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; check params
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(parlist) EQ 0 AND n_elements(addstruct) EQ 0 THEN BEGIN
      print,'Nothing to add!'
      return
  ENDIF 

  IF n_elements(parlist) EQ 0 THEN parlist = '' $
  ELSE parlist = strupcase(parlist)

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

  nf = n_elements(fnums)

  taglist = [ parlist, 'run', 'rerun', 'camcol', 'field', 'ID']
  taglist = taglist[ rem_dup(taglist) ]

  time=systime(1)

  FOR i=0, nf-1 DO BEGIN 
      
      field = fnums[i]
      outfile = newfiles[i]

      IF parlist[0] NE '' THEN BEGIN 
          struct=0
          read_tsobj, dir, struct, start=field, taglist=taglist, tsobjstr=tsobjstr
      ENDIF 
     
      corstruct = 0
      corrstruct = mrdfits(corrfiles[i], 1, hdr, /silent)

      ;; REMEMBER TO SWITCH BACK!!

      hdr = headfits(tsfiles[i])
      hdr1 = headfits(tsfiles[i], ext=1)
      FXHCLEAN,hdr              ;Mwrfits needs clean hdr
      SXADDPAR, hdr, 'SEEING_U', SXPAR(hdr1, 'SEEING_U')
      SXADDPAR, hdr, 'SEEING_G', SXPAR(hdr1, 'SEEING_G')
      SXADDPAR, hdr, 'SEEING_R', SXPAR(hdr1, 'SEEING_R')
      SXADDPAR, hdr, 'SEEING_I', SXPAR(hdr1, 'SEEING_I')
      SXADDPAR, hdr, 'SEEING_Z', SXPAR(hdr1, 'SEEING_Z')

      ncorr = n_elements(corrstruct)
      
      ;; check the tags
      IF i EQ 0 THEN BEGIN

          corrtags = tag_names(corrstruct)
          s = corrstruct[0]
          npar = 0
          nadd = 0

          IF parlist[0] NE '' THEN BEGIN ;; may have only given addstruct
              oldtags = tag_names(struct)
              npar = npar + n_elements(parlist)

              FOR kk=0, npar-1 DO BEGIN
                  ww=where(corrtags EQ parlist[kk], nww)
                  IF nww NE 0 THEN BEGIN 
                      print
                      print,'Tag ',parlist[kk],' is already in the structure'
                      return
                  ENDIF 
              ENDFOR 
          
              FOR kk=0, npar-1 DO BEGIN
                  wt = where(oldtags EQ parlist[kk])
                  nn = n_elements( struct[0].(wt[0]) )
                  s=create_struct(s, parlist[kk], struct[0].(wt[0]) )
              ENDFOR 
          ENDIF 

          ;; Now check addstruct
          IF n_elements(addstruct) NE 0 THEN BEGIN 
              addtags = tag_names(addstruct)
              nadd = n_elements(addtags)
              FOR kk=0, nadd-1 DO BEGIN 
                  ww=where(parlist EQ addtags[kk] OR corrtags EQ addtags[kk], nww)
                  IF nww NE 0 THEN BEGIN
                      print
                      print,'Add list tag ',addtags[kk],' is a duplicate'
                      return
                  ENDIF 
              ENDFOR
              
              s = create_struct(s, addstruct)
          ENDIF 

          IF npar EQ 0 AND nadd EQ 0 THEN BEGIN
              print,'Nothing to add!'
              return
          ENDIF 

          newtags = tag_names(s)

      ENDIF 

      IF parlist[0] NE '' THEN BEGIN 
          photo_match, struct.run, struct.rerun, struct.camcol, $
            struct.field, struct.id, $
            corrstruct.run, corrstruct.rerun, corrstruct.camcol, $
            corrstruct.field,corrstruct.id, $
            match, corrmatch

          ncorrmatch = n_elements(corrmatch)
          IF ncorrmatch NE ncorr THEN BEGIN
              print,'WHat!'
              return
          ENDIF 
      ENDIF 

      newstruct = 0
      newstruct = replicate(s, ncorr)

      copy_struct, corrstruct, newstruct

      FOR kk=0, npar-1 DO BEGIN 
          wold = where(oldtags EQ parlist[kk])
          wnew = where(newtags EQ parlist[kk])

          newstruct[corrmatch].(wnew[0]) = struct[match].(wold[0])
      ENDFOR 

      diff_list = compare_struct(corrstruct,newstruct)
      IF diff_list.ndiff NE 0 THEN BEGIN
          help,diff_list,/str
          return
      ENDIF 

;      help,newstruct,/str

      print,outfile
      mwrfits, newstruct, outfile, hdr, /create
;      return

      
  ENDFOR 
  print,'----------------'
  ptime,systime(1)-time

  return
END 
