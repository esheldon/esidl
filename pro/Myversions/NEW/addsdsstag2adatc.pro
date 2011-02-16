PRO addsdsstag2adatc, run, rerun, camcol, tags, overwrite=overwrite, ex_struct=ex_struct

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: addsdsstag2adatc, run, rerun, camcol, tags, overwrite=overwrite'
      return
  ENDIF 

  IF n_elements(ex_struct) NE 0 THEN BEGIN 
      do_ex_struct=1
      extags = tag_names(ex_struct)

      match, extags, tags, m1, m2
      IF m1[0] NE -1 THEN message,'ex_struct cannot contain tags common with input tags'
  ENDIF ELSE do_ex_struct=0

  front='adatc'

  IF NOT keyword_set(overwrite) THEN newfront = 'adatcadd' ELSE newfront=front

  fetch_dir, run, camcol, rerun, dir, corrdir=corrdir, /check

  IF corrdir NE '' THEN BEGIN 

      fetch_file_list, corrdir, files, fnums, front=front
  
      FOR i=0L, n_elements(fnums)-1 DO BEGIN 
          
          infile = files[i]
          outfile = repstr(infile, front, newfront)

          ;; can't use mrdfits3 with /deja_vu on at same time on
          ;; more than one type of file, so use regular mrdfits here
          tmp = mrdfits(infile, 1, hdr1)
          hdr0=headfits(infile, exten=0)
          fxhclean, hdr1

          IF ( size(tmp) )(0) EQ 0 THEN BEGIN 
              print,files[i] +' is empty'
              newstruct=tmp
          ENDIF ELSE BEGIN 
              addsdsstag, tmp, tags, newstruct
              IF do_ex_struct THEN BEGIN 
                  nn=n_elements(newstruct)
                  tex=replicate(ex_struct, nn)
                  combine_structs, temporary(newstruct), tex, tmpstruct
                  newstruct = temporary(tmpstruct)
              ENDIF 
          ENDELSE 

          print,'Output file: '+outfile
          mwrfits2, newstruct, outfile, hdr1, /create, /destroy, hdr0=hdr0

          tmp=0

      ENDFOR 

  ENDIF ELSE BEGIN 
      message,' no such run/rerun/camcol'
  ENDELSE 

END 
