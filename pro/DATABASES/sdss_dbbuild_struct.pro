PRO sdss_dbbuild_struct_open, database, create=create, noprompt=noprompt,$
                               silent=silent

  IF keyword_set(create) THEN BEGIN
      IF NOT keyword_set(noprompt) THEN BEGIN 
          ans = ' '
          read, ans, $
            prompt='Are you sure you want to initialize '+database+' (y/n)?'
      ENDIF ELSE ans = 'y'
      IF (ans EQ 'y') OR (ans EQ 'Y') THEN BEGIN 
          IF NOT keyword_set(silent) THEN $
            print,'Initializing database: '+database
          dbcreate, database, 1, 1
      ENDIF ELSE return
  ENDIF ELSE BEGIN 
      IF NOT keyword_set(silent) THEN print,'Appending database: '+database
  ENDELSE 
  dbopen,database,1

END 

PRO sdss_dbbuild_struct, struct, database, $
                          create=create, noindex=noindex, noprompt=noprompt,$
                          sort_photoid=sort_photoid, silent=silent

  IF  n_params() LT 2 THEN BEGIN 
      print,'-Syntax: sdss_dbbuild_struct, struct, database, /create, '
      print,'             /noindex, /noprompt, /sort_photoid'
      print
      print,'Assumes database was build with order of items=order of tags'
      print,'   plus photoid at the beginning'
      print,'Appends database.  set /create to Initialize database'
      print,'/sort will sort by photoid'
      return
  ENDIF 

  IF NOT keyword_set(silent) THEN BEGIN 
      print
      print,'Database: ',database
  ENDIF 

  !priv=2
  sdss_dbbuild_struct_open, database, $
    create=create, noprompt=noprompt, silent=silent

  phid = photoid(struct)

  IF keyword_set(sort_photoid) THEN BEGIN 
      IF NOT keyword_set(silent) THEN BEGIN 
          print
          print,'Sorting by photoid'
      ENDIF 
      s = sort(phid)

      struct = struct[s]
      phid = phid[s]
  ENDIF 

  ;; build the execution string
  tags = tag_names(struct)
  ntags = n_elements(tags)
  command = $
    'ess_dbbuild, phid'
  FOR i=0L, ntags-1 DO BEGIN 
      command = command + ', struct.'+tags[i]
  ENDFOR 
  command = command + ', noindex=noindex'

  IF NOT keyword_set(silent) THEN print,'Adding objects to database'

  IF NOT execute(command) THEN message,'Failure'

END 
