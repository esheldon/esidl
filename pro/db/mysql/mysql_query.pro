FUNCTION mysql_query, query, $
                      host=host, user=user, pwd=pwd, $
                      database=database, $
                      nrows=nrows, $
                      field_names=field_names, $
                      field_types=field_types, $
                      field_flags=field_flags, $
                      field_lengths=field_lengths, $
                      noconvert=noconvert, $
                      print_rows=print_rows, $
                      all=all, $
                      memory_check=memory_check, $
                      status=status

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: result=mysql_query(query, '
      print,'                            host=, user=, pwd=,'
      print,'                            database=, '
      print,'                            nrows=, '
      print,'                            field_names=, field_types=, '
      print,'                            field_flags=, field_lengths=,'
      print,'                            /noconvert, '
      print,'                            /print_rows, /all, '
      print,'                            /memory_check, '
      print,'                            status='
      return,-1
  ENDIF 

  ;; Send the query and return the rows of data
  rows = mysql_query_database(query,$
                              host=host, user=user, pwd=pwd, $
                              database=database, $
                              field_names=field_names, $
                              field_types=field_types, $
                              field_flags=field_flags, $
                              field_lengths=field_lengths, $
                              nrows=nrows, $
                              status=status)

  IF (status NE 0) THEN return,-1
  IF (nrows EQ 0) THEN return,-1
  IF keyword_set(noconvert) THEN return,rows

  IF keyword_set(memory_check) THEN help,/mem

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Build the structure.  This causes factor of ~2 overhead in memory
  ;; This is unavoidable unless we read line-by-line, but then we don't
  ;; know how many results are coming!
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  struct = mysql_getstruct(field_names, field_types, field_flags)
  struct = replicate(struct, nrows)
  IF keyword_set(memory_check) THEN help,/mem

  ;; Read into the structure
  reads, rows, struct
  IF keyword_set(memory_check) THEN help,/mem
  
  ;; Free the original data. Could do with temporary, but this
  ;; allows us to see the freed memory
  rows = 0
  IF keyword_set(memory_check) THEN help,/mem

  ;; Print some rows
  IF keyword_set(print_rows) THEN BEGIN 
      IF keyword_set(all) OR print_rows[0] GE nrows THEN BEGIN 
          print_struct, struct
      ENDIF ELSE BEGIN 
          nprint = print_rows[0] < nrows > 1
          print_struct, struct[0:nprint-1]
      ENDELSE 
  ENDIF 

  return,struct

END 
