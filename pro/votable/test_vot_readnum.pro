PRO test_vot_readnum

  dir = '~/tmp/VOlib_0.1/data/'

  ;; takes 1 min 32.2 sec with Chris'
  file = dir + 'tsObj-000756-6-44-0367.vot'
;  file = dir + 'votable.xml'

  tm = systime(1) 
  struct = vot_readnum(file)
  ptime,systime(1)-tm

END 
