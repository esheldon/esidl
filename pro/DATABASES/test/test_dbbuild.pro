PRO test_dbbuild, create=create, noprompt=noprompt

  !priv = 2

  ;database = '/net/cheops1/data4/sdss_database/adatc/adatc_test'
  database = '~/data/dbtest'

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

    n=1000
    stdef = {id:0L, ra:0d, dec:0d}
    cat=replicate(stdef, n)
    cat.id = lindgen(n)
    cat.ra = randomu(seed,n,/double)
    cat.dec = randomu(seed,n,/double)

    dbbuild, $
        cat.id, $
        cat.ra, $
        cat.dec

END 
