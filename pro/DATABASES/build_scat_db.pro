PRO build_scat_db, stripe, $
                   create=create, noindex=noindex, noprompt=noprompt

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: build_scat_db, stripe, clr, /create, /noindex, /noprompt'
      print
      print,'Appends database.  set /create to Initialize database'
      return
  ENDIF 

  IF stripe GE 9 AND stripe LE 15 THEN chunk=1
  IF stripe GE 27 AND stripe LE 37 THEN chunk=2
  IF stripe GE 42 AND stripe LE 86 THEN chunk=3

  cd,getenv("ZDBASE")

  clr = [1,2,3]
  clrstr = clrarr2string(clr)

  database = 'scat'+ntostr(chunk)
  print
  print,'Stripe: ',stripe
  print,'Color:  ',clrstr
  print,'Database: ',database
  print

  IF NOT keyword_set(noprompt) THEN BEGIN 
      ans = ' '
      read, ans, $
        prompt='Is this information correct (y/n)?'
      IF (ans NE 'y') AND (ans NE 'Y') THEN return
      print
  ENDIF 

  !priv=2
  
  IF keyword_set(create) THEN BEGIN
      IF NOT keyword_set(noprompt) THEN BEGIN 
          ans = ' '
          read, ans, $
            prompt='Are you sure you want to initialize '+database+' (y/n)?'
      ENDIF ELSE ans = 'y'
      IF (ans EQ 'y') OR (ans EQ 'Y') THEN BEGIN 
          print,'Initializing database: '+database
          dbcreate, database, 1, 1
      ENDIF ELSE return
  ENDIF ELSE print,'Appending database: '+database
  dbopen,database,1

  ;; Get the catalog
  IF n_elements(scat) EQ 0 THEN get_scat, stripe, clr, scat, /hirata
  count=n_elements(scat)

  photoz_select_struct, scat, tuseind

  useind = replicate(0b, count)
  IF tuseind[0] NE -1 THEN useind[tuseind] = 1b


  tl = tag_names(scat)
  w=where(tl EQ 'LAMBDA',nw)
  
  n=n_elements(scat)

  stripe_out = replicate(stripe, n)
  phid = photoid(scat)

  print
  print,'Adding objects to database'

  time = systime(1)

  ess_dbbuild, $
    phid, $
    stripe_out, $
    useind, $
    scat.run, $
    scat.rerun, $
    scat.camcol, $
    scat.field, $
    scat.id, $
    scat.e1, $
    scat.e2, $
    scat.e1_recorr, $
    scat.e2_recorr, $
    scat.e1e1err, $
    scat.e1e2err, $
    scat.e2e2err, $
    scat.upetro, $
    scat.gpetro, $
    scat.rpetro, $
    scat.ipetro, $
    scat.zpetro, $
    scat.grmodel, $
    scat.rimodel, $
    scat.rpetrorad, $
    scat.rsmear, $
    scat.clambda, $
    scat.ceta, $
    scat.objc_prob_psf, $
    scat.photoz_z, $
    scat.photoz_zerr, $
    scat.photoz_quality, $
    scat.photoz_abscounts, $    
    scat.nrunave, $
    scat.runcombine_flags, $
    scat.leafid, $
    noindex=noindex


  ptime,systime(1)-time

  return
END 
