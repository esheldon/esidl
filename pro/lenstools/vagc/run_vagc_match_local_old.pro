PRO run_vagc_match_local, overwrite=overwrite, lrg=lrg

  vagc_getstripes, stripes, nstripe

  IF keyword_set(lrg) THEN lrgstr = '_lrg' ELSE lrgstr = ''

  dir = vagc_lensinput_dir()

  taglist = $
    ['run','rerun','camcol','field','id','ra','dec']

  FOR i=0L, nstripe-1 DO BEGIN 

      sstr = stripe2string(stripes[i])
      infile = dir + 'stripe'+sstr+lrgstr+'_vagc.fit'
      outfile = dir + 'stripe'+sstr+lrgstr+'_vagc_matchlocal.fit'

      IF keyword_set(overwrite) OR (NOT fexist(outFile)) THEN BEGIN 
          print
          print,'Reading input file: ',inFile

          str = mrdfits(inFile,1)
          print
          print,'Will write to file: ',outfile
          
          vagc_match_local, str, taglist, match_struct
          
          add_tags, $
            temporary(str), ['match_rerun','match_id'], ['-1','-1L'], newstr
          
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
          
          print,'Writing file: ',outFile
          mwrfits, newstr, outFile, /create
          
      ENDIF 
  ENDFOR 

END 
