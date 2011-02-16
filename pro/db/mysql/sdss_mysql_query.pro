;; Same as mysql_query, but converts columns to sdss-style 5-arrays
;; when applicable.  The requested columns are sent through the columns
;; keyword rather than being part of the statement.  The statement is
;; the portion after the select columns from table part.  Table is also sent
;; through a keyword, with the default being tsObj.  Similarly, the default
;; database is sdss.  The default host not set here, but perhaps this is
;; a configured variable, or do we expect the user has set up .my.cnf as
;; for mysql_query?

FUNCTION sdss_mysql_query, statement, $
                           columns=columns, $
                           host=host, user=user, pwd=pwd, $
                           database=database, table=table, $
                           nrows=nrows, $
                           field_names=field_names, $
                           field_types=field_types, $
                           field_flags=field_flags, $
                           field_lengths=field_lengths, $
                           noconvert=noconvert, $
                           print_rows=print_rows, $
                           all=all, $
                           status=status

  status = 1

  IF n_elements(user) EQ 0 THEN user = 'sdss'
  IF n_elements(database) EQ 0 THEN database = 'sdss'
  IF n_elements(table) EQ 0 THEN table = 'tsObj'

  IF n_elements(columns) EQ 0 THEN BEGIN 

      ;; by default, we will read all columns.  We must query the database
      ;; to get the tags for this table

      tagsFromQuery = 1

      q = 'describe '+table
      colst = mysql_query(q, $
                          host=host, user=user, pwd=pwd, $
                          database=database, $
                          status=tstatus, nrows=tnrows)
      IF tstatus NE 0 THEN BEGIN  
          message,'Unable to retrieve column names',/inf
          return,-1
      ENDIF 

      tagObj = obj_new('mysql_taglist')
      tagObj->set, colst.field, 'MYSQL'

  ENDIF ELSE BEGIN

      tagsFromQuery = 0

      ;; Convert the sdss columns to mysql columns
      tagObj = obj_new('mysql_taglist')
      tagObj->set, columns, 'SDSS'
      tagObj->convert
      
  ENDELSE 

  taglist = strjoin(tagObj->tags(), ', ')
  query = 'SELECT '+taglist+' FROM '+table+' '+statement
  
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

  ;; convert the fields to sdss style, with arrays of 5
  tagObj->convert, sub

  ;; build the structure
  struct = $
    mysql_getstruct(tagObj->tags(), field_types[sub], field_flags[sub],$
                    nvalues = tagObj->nvalues())

  struct = replicate(struct, nrows)

  ;; Read into the structure
  reads, rows, struct

  ;; Print some rows
  IF keyword_set(print_rows) THEN BEGIN 
      IF keyword_set(all) OR print_rows[0] GE nrows THEN BEGIN 
          print_struct, struct
      ENDIF ELSE BEGIN 
          nprint = print_rows[0] < nrows > 1
          print_struct, struct[0:nprint-1]
      ENDELSE 
  ENDIF 

  status = 0
  return, struct

END 
