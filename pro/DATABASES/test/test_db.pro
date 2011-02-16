PRO test_db, struct

  dbopen, 'scat_gri'

  list = dbget('stripe',9)
  help,list
  dbclose
return
  list = 1+lindgen(10)
  dbprint, list, ['photoid','stripe','photoz_use','clambda','ceta']

  dbclose
  
return
  tm=systime(1)
  list = dbfind('id > 0')
  print,systime(1)-tm

  list = list[0:9]


  items = ['entry','id', 'ra', 'dec']
  items = ['entry','ra']
  ;;items = '*'
  tm=systime(1)
  dbprint,list, items
  print,systime(1)-tm

;  tm=systime(1)
;  dbext, list, items, id, ra, dec
;  print,systime(1)-tm

  dbclose

stop
END 
