PRO test_spec_dbbuild, cat, create=create, noprompt=noprompt

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: test_spec_dbbuild, cat, /create'
      return
  ENDIF 

  !priv = 2

  database = 'test_spec'

  IF keyword_set(create) THEN BEGIN
      IF NOT keyword_set(noprompt) THEN BEGIN 
          ans = ' '
          read, ans, $
            prompt='Are you sure you want to initialize '+database+' (y/n)?'
      ENDIF ELSE ans = 'y'
      IF (ans EQ 'y') OR (ans EQ 'Y') THEN BEGIN 
          print,'Initializing database: '
          dbcreate, database, 1, 1
      ENDIF ELSE return
  ENDIF ELSE print,'Appending database: '+database
  dbopen,database,1

  specid = sdss_specid(cat)

  s=sort(specid)
  specid = specid[s]
  cat = cat[s]

  depth = 10
  htmid = htm_index(cat.ra, cat.dec, depth)

  tm=systime(1)

  dbbuild, $
    specid, $
    htmid, $
    cat.plate, $
    cat.mjd, $
    cat.fiberid, $
    cat.z, $
    cat.z_err, $
    cat.zwarning, $
    cat.counts_model_pogson

  ptime,systime(1)-tm

  ;; No need to close, it is closed by dbbuild

END 
