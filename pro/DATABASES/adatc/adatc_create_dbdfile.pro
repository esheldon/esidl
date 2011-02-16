PRO adatc_create_dbdfile, run, camcol, rerun=rerun, database=database

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: adatc_create_dbdfile, run, camcol, rerun=, database='
      return
  ENDIF 

  database = adatc_dbfile(run,camcol, rerun=rerun, dir=outdir)
  database = outdir+database

  dbdfile = database + '.dbd'

  IF NOT fexist(outdir) THEN spawn,'mkdir '+outdir

  ;; max of 2 million
  maxentries = '2000000'

  ;; Read in one catalog to get tags
  cat = sdss_read('adatc',run, camcol, rerun=rerun)
  cat = cat[0]

  IF NOT tag_exist(cat,'htmindex') THEN BEGIN 
      cat = create_struct(cat, 'htmindex', 0L)
  ENDIF 
  
  db_struct2def, cat, names, typedefs

  print
  print,'Creating dbd file: ',dbdfile
  openw, lun, dbdfile, /get_lun

  printf,lun,'#title'
  printf,lun,'Adatc catalog database for run '+$
    ntostr(run)+' camcol '+ntostr(camcol)
  printf,lun
  printf,lun,'#maxentries'
  printf,lun, maxentries
  printf,lun
  printf,lun,'#items'
  
  ;; First some extra items
  printf,lun,'photoid                U*8     SDSSid'
  maxlen = max(strlen(names))
  FOR i=0L, n_elements(names)-1 DO BEGIN 

      tn = names[i]
      tnl = strlen(tn)
      IF tnl LT maxlen THEN BEGIN 
          tn = tn + strjoin(replicate(' ', maxlen-tnl))
      ENDIF 
      printf,lun,tn,typedefs[i],format='(a,a7)'

  ENDFOR 

  printf,lun
  printf,lun,'#formats'
  printf,lun
  printf,lun,'#index'
  printf,lun,'photoid    sorted'
  ;printf,lun,'spec_entry index'
  printf,lun,'htmindex     sort'
  printf,lun,'field        sort'


  free_lun, lun

END 
