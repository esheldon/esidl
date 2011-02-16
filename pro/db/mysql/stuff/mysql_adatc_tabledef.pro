PRO mysql_adatc_tabledef, sqlfile

  IF n_elements(sqlfile) EQ 0 THEN BEGIN 
      lun = -1
  ENDIF ELSE BEGIN 
      openw,lun,sqlfile,/get_lun
  ENDELSE 
  
  struct = make_corrected_struct()
  coldefs = mysql_struct2coldefs(struct, /sdss, outtags=outtags)



  w=where(outtags EQ 'rerun')
  coldefs[w[0]] = 'rerun SMALLINT NOT NULL'
  w=where(outtags EQ 'camcol')
  coldefs[w[0]] = 'camcol SMALLINT NOT NULL'
  w=where(outtags EQ 'field')
  coldefs[w[0]] = 'field SMALLINT NOT NULL'

  IF 1 THEN BEGIN 
      coldefs = [coldefs, $
                 'PRIMARY KEY (photoid)', $
                 'UNIQUE INDEX ind_rrcfi (run,rerun,camcol,field,id)', $
                 $
                 'INDEX ind_value_flags (value_flags)', $
                 'INDEX ind_corrselect_flags (corrselect_flags)', $
                 $
                 'INDEX ind_prob_gal (objc_prob_gal)', $
                 $
                 'INDEX ind_cmodel_r (cmodel_counts_r)', $
                 'INDEX ind_cmodel_ext_r (cmodel_counts_ext_r)', $
                 $
                 'INDEX ind_m_r_r (m_r_r)', $
                 'INDEX ind_m_r_h_r (m_r_h_r)', $
                 'INDEX ind_seeing_r (seeing_r)', $
                 'INDEX ind_htm_index (htm_index)']

      ncoldefs = n_elements(coldefs)

      printf,lun,'use sdss;'
      printf,lun,'CREATE TABLE adatc'
      printf,lun,'('
      FOR i=0L, ncoldefs-2 DO BEGIN 
          printf,lun,coldefs[i]+', '
      ENDFOR 
      printf,lun,coldefs[i]
      printf,lun,') MAX_ROWS = 500000000;'

      IF n_elements(sqlfile) NE 0 THEN free_lun, lun

  ENDIF ELSE BEGIN 
      coldefs = [coldefs, $
                 'PRIMARY KEY (photoid)', $
                 'UNIQUE INDEX ind_rrcfi (run,rerun,camcol,field,id)', $
                 'INDEX ind_cmodel_r (cmodel_counts_r)', $
                 'INDEX ind_cmodel_ext_r (cmodel_counts_ext_r)', $
                 'INDEX ind_m_r_r (m_r_r)', $
                 'INDEX ind_m_r_h_r (m_r_h_r)', $
                 'INDEX ind_seeing_r (seeing_r)', $
                 'INDEX ind_prob_gal (objc_prob_gal)', $
                 'INDEX ind_htm_index (htm_index)']

      ncoldefs = n_elements(coldefs)

      printf,lun,'use sdss;'
      printf,lun,'CREATE TABLE adatc'
      printf,lun,'('
      FOR i=0L, ncoldefs-2 DO BEGIN 
          printf,lun,coldefs[i]+', '
      ENDFOR 
      printf,lun,coldefs[i]
      printf,lun,') MAX_ROWS = 500000000;'

      IF n_elements(sqlfile) NE 0 THEN free_lun, lun

  ENDELSE 


END 
