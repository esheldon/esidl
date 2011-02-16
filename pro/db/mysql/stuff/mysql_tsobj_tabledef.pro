PRO mysql_tsobj_tabledef, sqlfile

  IF n_elements(sqlfile) EQ 0 THEN BEGIN 
      lun = -1
  ENDIF ELSE BEGIN 
      openw,lun,sqlfile,/get_lun
  ENDELSE 
  phototags = mysql_tsobj_stuff_taglist(struct=struct)

  coldefs = mysql_struct2coldefs(struct, /sdss, outtags=outtags)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; No need for int32 on these columns
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  w=where(outtags EQ 'rerun')
  coldefs[w[0]] = 'rerun SMALLINT NOT NULL'
  w=where(outtags EQ 'camcol')
  coldefs[w[0]] = 'camcol SMALLINT NOT NULL'
  w=where(outtags EQ 'field')
  coldefs[w[0]] = 'field SMALLINT NOT NULL'
  w=where(outtags EQ 'parent')
  coldefs[w[0]] = 'parent SMALLINT NOT NULL'

; Crowded fields in segue runs? may need int32
; max id in database is 3517
;  w=where(outtags EQ 'id')
;  coldefs[w[0]] = 'id SMALLINT NOT NULL'

  w=where(outtags EQ 'nchild')
  coldefs[w[0]] = 'nchild SMALLINT NOT NULL'
  w=where(outtags EQ 'objc_type')
  coldefs[w[0]] = 'objc_type SMALLINT NOT NULL'

  w=where(outtags EQ 'type_u')
  coldefs[w[0]] = 'type_u SMALLINT NOT NULL'
  w=where(outtags EQ 'type_g')
  coldefs[w[0]] = 'type_g SMALLINT NOT NULL'
  w=where(outtags EQ 'type_r')
  coldefs[w[0]] = 'type_r SMALLINT NOT NULL'
  w=where(outtags EQ 'type_i')
  coldefs[w[0]] = 'type_i SMALLINT NOT NULL'
  w=where(outtags EQ 'type_z')
  coldefs[w[0]] = 'type_z SMALLINT NOT NULL'

  w=where(outtags EQ 'nprof_u')
  coldefs[w[0]] = 'nprof_u SMALLINT NOT NULL'
  w=where(outtags EQ 'nprof_g')
  coldefs[w[0]] = 'nprof_g SMALLINT NOT NULL'
  w=where(outtags EQ 'nprof_r')
  coldefs[w[0]] = 'nprof_r SMALLINT NOT NULL'
  w=where(outtags EQ 'nprof_i')
  coldefs[w[0]] = 'nprof_i SMALLINT NOT NULL'
  w=where(outtags EQ 'nprof_z')
  coldefs[w[0]] = 'nprof_z SMALLINT NOT NULL'

  ;; Fix dec reserved word problem
  w=where(outtags EQ 'dec')
  coldefs[w[0]] = 'decl DOUBLE NOT NULL'

  coldefs = ['photoid BIGINT UNSIGNED NOT NULL', $
             coldefs, $
             'htm_index BIGINT UNSIGNED NOT NULL', $
             'PRIMARY KEY (photoid)']

  ncoldefs = n_elements(coldefs)

  printf,lun,'use sdss;'
  printf,lun,'CREATE TABLE tsObj'
  printf,lun,'('
  FOR i=0L, ncoldefs-2 DO BEGIN 
      printf,lun,coldefs[i]+', '
  ENDFOR 
  printf,lun,coldefs[i]
  printf,lun,') MAX_ROWS = 500000000;'

  IF n_elements(sqlfile) NE 0 THEN free_lun, lun
END 
