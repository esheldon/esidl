PRO combine_zlensum_files, files, struct, count=count, columns=columns, silent=silent

  ;; adds up the zindex's
  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: combine_zlensum_files, files, struct'
      return
  ENDIF 

  nfiles = n_elements(files)

  ;; first go through, read the headers and find how many objects
  ;; are in each.  Then allocate a big struct and copy it in.

  IF NOT keyword_set(silent) THEN BEGIN
      print
      print,'Reading headers'
  ENDIF 
  ntotal = 0L
  numlist = lonarr(nfiles)
  FOR i=0L, nfiles-1 DO BEGIN 

      hdr = headfits(files[i], ext=1)
      numlist[i] = sxpar(hdr,'naxis2')
      ntotal = ntotal + numlist[i]

  ENDFOR 

  IF NOT keyword_set(silent) THEN BEGIN
      print
      print,'Total number of objects: '+ntostr(ntotal)
  ENDIF 

  rnd = ntostr(long(randomu(seed)*10000L))
  typ = 'tmp'+rnd+'_'+ntostr(long(systime(1)))
  beg =0L
  FOR i=0L, nfiles-1 DO BEGIN 

      IF numlist[i] NE 0 THEN BEGIN 
          print,'Reading File: '+ntostr(files[i])
          t=mrdfits(files[i], 1, hdr, silent=silent, columns=columns)
          IF i EQ 0 THEN BEGIN 

              struct_tags = tag_names(t[0])
              struct_ind  = lindgen(n_elements(struct_tags))
              
              ;; "zero" the struct for those tags that
              ;; will not be copied.
              struct = create_struct(t[0], name=typ)
              zero_struct,struct
              
              ;; simple copy on first one
              struct=replicate(struct, ntotal)
              struct[beg:beg+numlist[i]-1] = t
              beg = beg+numlist[i]
              
          ENDIF ELSE BEGIN 

              max_zindex = max(struct.zindex)          
              t.zindex = t.zindex + (max_zindex + 1)

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; must match up the tags for possibly different struct
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              newtags = tag_names(t[0])
              match, struct_tags, newtags, mstr, mnew
              IF mstr[0] EQ -1 THEN $
                message,'No compatible tags from file: '+files[i]
              
              nmatch = n_elements(mnew)
              ;; loop over tags
              FOR j=0L, nmatch-1 DO BEGIN 
                  si = mstr[j]
                  sn = mnew[j]
                  struct[beg:beg+numlist[i]-1].(si) = t.(sn)
              ENDFOR 
              
              beg = beg+numlist[i]

          ENDELSE 
          
          t = 0
          
      ENDIF ELSE BEGIN  
          print,'File is empty: '+ntostr(files[i])
      ENDELSE 
  ENDFOR 

  count = ntotal

return


END 
