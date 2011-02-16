PRO build_spectra_db, stripe_in, create=create, noprompt=noprompt

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: build_spectra_db, stripe_in, create=create, noprompt=noprompt'
      print
      print,'Appends database.  set /create to Initialize database'
      print,'You must run this in the $ZDBASE directory!!!'
      return
  ENDIF 

  ;; You must run this in the $ZDBASE directory!!!
  time = systime(1)
  database = 'spectra'
  IF NOT keyword_set(noprompt) THEN noprompt=0
  print
  print,'Stripe: ',stripe_in
  print,'Database: ',database
  print
  IF NOT noprompt THEN BEGIN 
      ans = ' '
      read, ans, $
        prompt='Is this information correct (y/n)?'
      IF (ans NE 'y') AND (ans NE 'Y') THEN return
      print
  ENDIF 

  get_spectra_lcat, stripe_in, str

  !priv=2
  print
  IF keyword_set(create) THEN BEGIN
      IF NOT noprompt THEN BEGIN 
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
  
  n=n_elements(str)
  stripe = replicate(stripe_in, n)



  dbbuild, $
    stripe,$
    str.run,$
    str.rerun,$
    str.camcol,$
    str.field,$
    str.id,$
    str.parent,$
    str.nchild,$
    str.type,$
    str.objc_rowc,$
    str.objc_colc,$
    str.rowc,$
    str.colc,$
    str.counts_model,$
    str.petrocounts,$
    str.petrorad,$
    str.status,$
    str.ra,$
    str.dec,$
    str.primtarget,$
    str.sectarget,$
    str.seeing,$
    str.e1,$
    str.e2,$
    str.momerr,$
    str.classification,$
    str.class,$
    str.z1d,$
    str.absmag,$
    str.lum

  ptime,systime(1)-time

  return
END 
