PRO test_imaging_dbbuild, cat, create=create, noprompt=noprompt

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: test_imaging_dbbuild, cat, /create'
      return
  ENDIF 

  ;; First, match up to the spec database
  ncat = n_elements(cat)
  w=where(cat.plate GT 0)

  specid = sdss_specid(cat[w])
  spec_entry = lonarr(ncat)

  dbopen,'test_spec'
  tspec_entry = dbmatch('specid', specid)
  dbclose

  spec_entry[w] = tspec_entry
  
  !priv = 2

  database = 'test_imaging'

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

  id = photoid(cat)

  s=sort(id)

  id = id[s]
  spec_entry = spec_entry[s]
  cat = cat[s]

  depth = 10
  htmid = htm_index(cat.ra, cat.dec, depth)

  tm=systime(1)

  ess_dbbuild, $
    id, $
    htmid, $
    $
    cat.run, $
    cat.rerun, $
    cat.field, $
    cat.camcol, $
    cat.id, $
    $
    cat.primtarget, $
    fix(cat.objc_type), $
    $
    cat.counts_model_pogson, $
    cat.counts_modelerr_pogson, $
    cat.objc_prob_psf, $
    $
    cat.photoz_z, $
    cat.photoz_zerr1, $
    cat.photoz_zerr2, $
    cat.photoz_zerr3, $
    fix(cat.photoz_zwarning), $
    $
    cat.ra, $
    cat.dec, $
    $
    spec_entry

  ptime, systime(1)-tm

  ;; No need to close, it is closed by dbbuild

return
  tm=systime(1)

  dbbuild, $
    id, $
    cat.run, $
    cat.rerun, $
    cat.field, $
    cat.camcol, $
    cat.id, $
    cat.ra, $
    cat.dec, $
    cat.photoz_z, $
    cat.photoz_zerr1, $
    cat.photoz_zerr2, $
    cat.photoz_zerr3, $
    cat.fiberid, $
    spec_entry

  ptime, systime(1)-tm

END 
