PRO run_vagc_match_local, overwrite=overwrite

  ;; Don't match lss files. Will store all in database and just
  ;; copy in the lss info.

  ;; In fact, the lss info should just be written out.

  dir = sdssidl_config('spec_dir')+'blanton/gal_collated/byrun/'

  vagc_files = file_search(dir+'run??????-vagc_collate.fits')
  lss_files = file_search(dir+'run??????-lss_sample14full0_collate.fits')

  files = [vagc_files, lss_files]
  nfiles = n_elements(files)

  FOR i=0L, nfiles-1 DO BEGIN 

      inFile = files[i]
      outFile = repstr(inFile, '.fits', '_matchlocal.fits')

;      IF 0 THEN BEGIN 
      IF keyword_set(overwrite) OR (NOT fexist(outFile)) THEN BEGIN 
          print
          print,'Reading input file: ',inFile

          str = mrdfits(inFile,1)
          print
          print,'Will write to file: ',outfile
          
          match_struct = vagc_match_local(str, status=status)
          
          add_tags, $
            temporary(str), ['match_rerun','match_id'], ['-1','-1L'], newstr
          
          ;; Found no matches if status ne 0
          IF status EQ 0 THEN BEGIN 
              ;; Do a check
              w=where( (newstr[match_struct.vagc_index].run NE $
                        match_struct.run) OR $
                       (newstr[match_struct.vagc_index].camcol NE $
                        match_struct.camcol) OR $
                       (newstr[match_struct.vagc_index].field NE $
                        match_struct.field), nw)
              
              IF nw NE 0 THEN BEGIN 
                  message,'Some do not match!!'
              ENDIF 
              
              newstr[match_struct.vagc_index].match_rerun = match_struct.rerun
              newstr[match_struct.vagc_index].match_id = match_struct.id
          ENDIF  
          
          print,'Writing file: ',outFile
          mwrfits, newstr, outFile, /create
          
      ENDIF 
  ENDFOR 

END 
