PRO build_lcat_db, lcat, stripe_in, clr_in, create=create, noprompt=noprompt

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: build_lcat_db, lcat, stripe, clr, create=create, noprompt=noprompt'
      print
      print,'Appends database.  set /create to Initialize database'
      print,'You must run this in the $ZDBASE directory!!!'
      return
  ENDIF 

  ;; You must run this in the $ZDBASE directory!!!

  time = systime(1)
  colors = ['u','g','r','i','z']
  IF NOT keyword_set(noprompt) THEN noprompt=0

  IF (clr_in[0] NE 1) AND (clr_in[0] NE 2) AND (clr_in[0] NE 3) THEN BEGIN 
      print,'clr must be 1,2, or 3'
      return
  ENDIF 

  CASE clr_in[0] OF 
      1: clr=clr_in[0]
      2: clr=clr_in[0]
      3: clr=clr_in[0]
      ELSE: BEGIN 
          print,'clr must be 1,2, or 3'
      END 
  ENDCASE 
  database = 'lcat_'+colors[clr]
  print
  print,'Stripe: ',stripe_in
  print,'Color:  ',clr,'  (',colors[clr],')'
  print,'Database: ',database
  print
  IF NOT noprompt THEN BEGIN 
      ans = ' '
      read, ans, $
        prompt='Is this information correct (y/n)?'
      IF (ans NE 'y') AND (ans NE 'Y') THEN return
      print
  ENDIF 

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

  tl = tag_names(lcat)
  w=where(tl EQ 'LAMBDA',nw)
  
  print
  IF nw EQ 0 THEN BEGIN 
      w=where(tl EQ 'RA', nw)
      IF nw EQ 0 THEN BEGIN 
          print,'Catalog must contain either RA,DEC or LAMBDA, ETA'
          return
      ENDIF ELSE BEGIN 
          print,'Converting ra and dec to lambda, eta'
          print
          eq2survey,lcat.ra,lcat.dec,lambda,eta
      ENDELSE 
  ENDIF ELSE BEGIN 
      eta = lcat.eta
      lambda = lcat.lambda
  ENDELSE 

  s = sort(lambda)
  n=n_elements(lcat)

  stripe = replicate(stripe_in[0], n)

  print,'Loading Database: '+database
  print

  stripe = replicate(stripe_in[0], n)

  dbbuild, $
    stripe, $
    lcat.petrocounts[2], $
    lcat.gr, $
    lcat.ri, $
    lcat.petror90, $
    lcat.e1, $
    lcat.e2, $
    lcat.momerr, $
    lcat.r, $
    lcat.classification, $
    lcat.photoz, $
    lcat.ra/15., $
    lcat.dec, $
    lambda, $
    eta

  ptime,systime(1)-time

  return
END 
