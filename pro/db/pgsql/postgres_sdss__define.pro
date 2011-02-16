FUNCTION postgres_sdss::init, connect_info=connect_info
  
  IF NOT self->postgres::init(connect_info=connect_info) THEN BEGIN 
      message,'Failed to initialize postgres',/inf
      return,0
  ENDIF 
  return,1

END 

; Execute the input set of queries.

FUNCTION postgres_sdss::send_queries, table, where_statements, $
                      columns=columns, $
                      clauses=clauses, $
                      connect_info=connect_info, $
                      slow=slow, $
                      status=status
  

  nQuery = n_elements(where_statements)
  ntable = n_elements(table)

  IF ntable NE 1 OR nquery EQ 0 THEN BEGIN 
      print,'-Syntax: struct = ps->send_queries(table, where_statements, '
      print,'                                   columns=, clauses=, '
      print,'                                   connect_info=, '
      print,'                                   /slow, status=)'
      print
      message,'Halting'
  ENDIF 

  IF n_elements(columns) NE 0 THEN BEGIN 
      cols = strjoin( columns )
  ENDIF ELSE BEGIN 
      cols = '*'
  ENDELSE 

  IF n_elements(clauses) EQ 0 THEN clauses = ''

  printlen = 124

  ;; Simple for a single query
  IF nQuery EQ 1 THEN BEGIN 

      query = $
        "SELECT "+cols+" FROM "+table+" WHERE "+where_statements+" "+clauses[0]
      print,strmid(query, 0, printlen)

      struct = pgsql_query(query, connect_info=connect_info, status=status)
      return,struct
  ENDIF 

  ;; Multiple Queries sent
  IF NOT keyword_set(slow) THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Fast: we gather all the objects in pointers and then
      ;; copy them into an output structure.  This is faster than the 
      ;; standard method below, but can use more memory for many separate
      ;; queries: a factor of two in the limit of an infinite number of
      ;; queries.   Its only slightly faster for large tables, however, when
      ;; grabbing the data takes much longer than the count(*) statement
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      queries = $
        "SELECT "+cols+" FROM "+table+" WHERE "+where_statements+" "+clauses[0]

      ptrlist = ptrarr(nQuery)
      nrows = 0ULL
      
      FOR i=0L, nQuery-1 DO BEGIN 

          print,strmid(queries[i], 0, printlen)
          st = pgsql_query(queries[i], nrows=tnrows, connect_info=connect_info, $
                           status=status)

          IF tnrows GT 0 THEN BEGIN 
              nrows = nrows + tnrows
              ptrlist[i] = ptr_new(st, /no_copy)
          ENDIF 

      ENDFOR 
      IF nrows GT 0 THEN BEGIN 
          struct = combine_ptrlist(ptrlist)
      ENDIF ELSE BEGIN 
          struct = -1
      ENDELSE 
      return, struct

  ENDIF ELSE BEGIN

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Standard: Count the queries results first by running a 
      ;; SELECT count(*) statement first, then create output and
      ;; copy in as we re-run queries.  For many large queries this 
      ;; could save lots of memory, but could also be a lot slower. 
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      numlist = ulon64arr(nQuery)
      nrows = 0ULL

      num_queries = $
        "SELECT count(*) from "+table+" WHERE "+where_statements+" "+clauses[0]

      FOR i=0L, nQuery-1 DO BEGIN 


          print, strmid(num_queries[i], 0, printlen)
          numstruct = pgsql_query(num_queries[i], nrows=cnrows, connect_info=connect_info, $
                                  status=status)

          IF cnrows NE 1 THEN BEGIN 
              tnrows = numstruct.count
              IF strnumber(tnrows) THEN BEGIN 
                  tnrows = ulong64(tnrows)
                  nrows = nrows + tnrows
                  numlist[i] = tnrows
              ENDIF 
          ENDIF ELSE BEGIN 
              return, -1
          ENDELSE 
      ENDFOR 

      ;; Any found?
      w=where(numlist GT 0, nw)
      IF nw NE 0 THEN BEGIN 

          ;; Now get a single structure
          struct_query = "SELECT "+cols+" FROM "+table+$
            " WHERE "+where_statements[w[0]]+" LIMIT 1"

          print,'getting struct'
          print,strmid(struct_query, 0, printlen)
          st = pgsql_query(struct_query, connect_info=connect_info, $
                           status=status)
          
          ;; Now all
          queries = $
            "SELECT "+cols+" FROM "+table+$
            " WHERE "+where_statements+" "+clauses[0]

          beg = 0ULL
          struct = replicate(st, nrows)
          FOR i=0ULL, nQuery-1 DO BEGIN 

              IF numlist[i] NE 0 THEN BEGIN 

                  print,strmid(queries[i], 0, printlen)
                  st = pgsql_query(queries[i], connect_info=connect_info, $
                                   status=status)
                  
                  struct[beg:beg+numlist[i]-1] = temporary(st)
                  beg=beg+numlist[i]
                  
              ENDIF 
          ENDFOR 
          
      ENDIF ELSE BEGIN 
          ;; None found
          return, -1
      ENDELSE 
      
  ENDELSE ;; Slower, memory efficient way


END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Individual ids requested
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION postgres_sdss::read_photoids, table, photoids, $
                      columns=columns, $
                      clauses=clauses, $
                      connect_info=connect_info, $
                      slow=slow

  IF n_elements(photoids) EQ 0 OR n_elements(table) EQ 0 THEN BEGIN 
      print,'-Syntax: struct = ps->read_photoids(table, photoids, '
      print,'                                    columns=, clauses=, '
      print,'                                    connect_info=, '
      print,'                                    /slow, status=)'
      print
      message,'Halting'
  ENDIF 

  ;; increased max_stack_depth to 16384 to handle large queries.
  maxnum = 10000

  tphotoids = photoids[rem_dup(photoids)]

  nid = n_elements(tphotoids)
  IF nid LT maxnum THEN BEGIN 

      photoidstr = ntostr(tphotoids)
      photoidstr = '('+strjoin(photoidstr, ',')+')'

      where_statement = 'photoid IN '+photoidstr
  ENDIF ELSE BEGIN 
      ndiv = nid/maxnum
      nmod = nid MOD maxnum
      nstatement = ndiv + (nmod NE 0)
      photoidstr = strarr(nstatement)

      beg = 0L
      FOR i=0L, nstatement-1 DO BEGIN 
          IF i EQ nstatement-1 THEN BEGIN 
              tid = tphotoids[beg:nid-1]
          ENDIF ELSE BEGIN 
              tid = tphotoids[beg:beg+maxnum-1]
          ENDELSE 
          tphotoidstr = ntostr(tid)
          photoidstr[i] = '('+strjoin(tphotoidstr, ',')+')'    

          tid = 0
          tphotoidstr = 0
          beg = beg + maxnum
      ENDFOR 

      ;; Now an array of statements
      where_statement = 'photoid IN '+photoidstr
  ENDELSE 
  struct = self->send_queries(table,where_statement, $
                              columns=columns, $
                              clauses=clauses, $
                              connect_info=connect_info, $
                              slow=slow, $
                              status=status)
  return,struct

END 

FUNCTION postgres_sdss::read_ids, table, runs, camcols, fields, ids, $
                      reruns=reruns_in, minrerun=minrerun, $
                      columns=columns, $
                      clauses=clauses, $
                      connect_info=connect_info, $
                      slow=slow

  
  IF n_params() LT 5 THEN BEGIN 
      print,'-Syntax: struct = ps->read_ids(table, runs, camcols, fields, ids, '
      print,'                               reruns=, '
      print,'                               columns=, clauses=, '
      print,'                               connect_info=, '
      print,'                               /slow, status=)'
      print
      message,'Halting'
  ENDIF 



  ;; always returns same length as runs
  reruns = $
    self->_get_reruns(runs, reruns=reruns_in, minrerun=minrerun)

  photoids = sdss_photoid(runs,reruns,camcols,fields,ids)
  IF photoids[0] EQ -1 THEN message,'Halting'

  return,self->read_photoids(photoids, slow=slow)

END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Field-by-field
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION postgres_sdss::read_fields, table, runs, camcols, fields, $
                      reruns=reruns_in, minrerun=minrerun, $
                      columns=columns, $
                      clauses=clauses, $
                      connect_info=connect_info, $
                      slow=slow

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: struct = ps->read_fields(table, runs, camcols, fields, '
      print,'                                  reruns=, '
      print,'                                  columns=, clauses=, '
      print,'                                  connect_info=, '
      print,'                                  /slow, status=)'
      print
      message,'Halting'
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; For non-consecutive fields, faster way gives a 30% improvement for 100 fields
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;; always returns same length as runs
  reruns = $
    self->_get_reruns(runs, reruns=reruns_in, minrerun=minrerun)

  nruns = n_elements(runs)
  ncamcols = n_elements(camcols)
  nfields = n_elements(fields)

  Ufields = fields[rem_dup(fields)]
  nUfields = n_elements(Ufields)

  IF nruns EQ 1 AND ncamcols EQ 1 AND nUfields EQ 1 THEN BEGIN 

      where_statements = $
        "run = "+ntostr(runs)+" AND " + $
        "rerun = "+ntostr(reruns) + " AND " + $
        "camcol = "+ntostr(camcols)+" AND "+ $
        "field = "+ntostr(Ufields)

  ENDIF ELSE IF nruns EQ 1 AND ncamcols EQ 1 THEN BEGIN 
      
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; single run, camcol, but more than one field
      ;; Detect the consecutive fields, which can speed things
      ;; up a lot
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;; Extract the consecutive ones
      get_consecutive, uFields, first, last
      nConsec = n_elements(first)

      ;; Generate where statements
      where_statements = strarr(nConsec)
      where_statements = $
        "run = "+ntostr(runs)+" AND " + $
        "rerun = "+ntostr(reruns)+" AND "+$
        "camcol = "+ntostr(camcols)+" AND "+ $
        "field BETWEEN "+$
        ntostr(uFields[first]) + " AND " + ntostr(uFields[last])

  ENDIF ELSE BEGIN

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Multiple runs,camcols.  User must input exactly same number of 
      ;; camcols and fields as runs in this situation
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF (ncamcols NE nruns) OR (nfields NE nruns) THEN BEGIN 
          message,'Number of camcols and fields must equal number '+$
            'of runs for multi-run, multi-camcol statements',/inf

          return,-1
      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Faster way, uses speed up for consecutive fields
      ;; For random fields is a 20-30% speedup.  Truly consecutive fields,
      ;; its a huge speedup
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ;; Make tree structure of run,rerun,camcol,field
      idstruct = sdss_histid(runs, reruns, camcols, fields)
      n_unique = idStruct.Nleaves
      ptrlist = ptrarr(n_unique)
      ptrIndex = 0L
      
      pruns = idStruct.runs
      FOR ri=0L, idStruct.nruns-1 DO BEGIN 
          runStr = ntostr( (*pruns)[ri].run )
          preruns = (*pruns)[ri].reruns
          FOR rri=0L, (*pruns)[ri].nreruns-1 DO BEGIN 
              rerunStr = ntostr( (*preruns)[rri].rerun )
              pcamcols = (*preruns)[rri].camcols
              FOR ci=0L, (*preruns)[rri].ncamcols-1 DO BEGIN 
                  camcolStr = ntostr( (*pcamcols)[ci].camcol )
                  fieldStructs = *(*pcamcols)[ci].fields
                  
                  rmd = rem_dup(fieldStructs.field)
                  Ufields = fieldStructs[rmd].field
                  nUfields = n_elements(Ufields)
                  
                  ;; Extract the consecutive ones
                  get_consecutive, uFields, first, last
                  nConsec = n_elements(first)
                  
                  ;; Generate where statements
                  twhere_statements = strarr(nConsec)
                  twhere_statements = $
                    "run = "+runStr+" AND " + $
                    "rerun = "+rerunStr+" AND "+$
                    "camcol = "+camcolStr+" AND "+ $
                    "field BETWEEN "+$
                    ntostr(uFields[first])+" AND "+ntostr(uFields[last])

                  ;; Now copy these into our pointerlist
                  ptrlist[ptrIndex] = ptr_new(twhere_statements, /no_copy)
                  ptrIndex = ptrIndex + 1
              ENDFOR ;; camcols
          ENDFOR ;; reruns
      ENDFOR                    ;runs
      
      ;; Now combine the lists
      where_statements = combine_ptrlist(ptrlist)
      sdss_histid_destroy, idstruct

  ENDELSE ;; multiple runs....
  
  struct = self->send_queries(table,where_statements, $
                              columns=columns, $
                              clauses=clauses, $
                              connect_info=connect_info, $
                              slow=slow, $
                              status=status)
  return,struct


END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Entire camcols
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION postgres_sdss::read_camcols, table, runs, camcols, $
                      reruns=reruns_in, minrerun=minrerun, $
                      columns=columns, $
                      clauses=clauses, $
                      connect_info=connect_info, $
                      slow=slow

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: struct = ps->read_camcols(table, runs, camcols, '
      print,'                                   reruns=, '
      print,'                                   columns=, clauses=, '
      print,'                                   connect_info=, '
      print,'                                   /slow, status=)'
      print
      message,'Halting'
  ENDIF 



  ;; always returns same length as runs
  reruns = $
    self->_get_reruns(runs, reruns=reruns_in, minrerun=minrerun)


  nruns = n_elements(runs)

  IF nruns EQ 1 THEN BEGIN 
      ;; Single run, perhaps multiple camcols
      where_statements = $
        "run = "+ntostr(runs)+" AND "+$
        "rerun = "+ntostr(reruns)+" AND "+$
        "camcol = "+ntostr(camcols)

  ENDIF ELSE BEGIN 

      ;; Multiple runs.  User must input exactly same number of camcols as
      ;; runs in this situation

      ncamcols = n_elements(camcols)
      IF ncamcols NE nruns THEN BEGIN 
          message,'Number of camcols must equal number of runs for '+$
            'multi-run statements',/inf

          return,-1
      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Not much speed advantage to optimizing this, so just split into
      ;; multiple unique camcol reads. Also, will use less memory to read
      ;; the individual camcols, since a copy is required. Note: doing a 
      ;; camcol IN (...) is actually SLOWER than doing the individual queries,
      ;; but still way faster than an OR on the camcols.
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      where_statements = $
        "run = "+ntostr(runs)+" AND "+$
        "rerun = "+ntostr(reruns)+" AND "+$
        "camcol = "+ntostr(camcols)
      where_statements = where_statements[ rem_dup(where_statements) ]

  ENDELSE 
  struct = self->send_queries(table,where_statements, $
                              columns=columns, $
                              clauses=clauses, $
                              connect_info=connect_info, $
                              slow=slow, $
                              status=status)
  return,struct


END 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Entire runs are requested.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION postgres_sdss::read_runs, table, runs, $
                      reruns=reruns_in, minrerun=minrerun, $
                      columns=columns, $
                      clauses=clauses, $
                      connect_info=connect_info, $
                      slow=slow

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: struct = ps->read_runs(table, runs, camcols, '
      print,'                                reruns=, '
      print,'                                columns=, clauses=, '
      print,'                                connect_info=, '
      print,'                                /slow, status=)'
      print
      message,'Halting'
  ENDIF 


  status=1

  ;; always returns same length as runs
  reruns = $
    self->_get_reruns(runs, reruns=reruns_in, minrerun=minrerun)

  nruns = n_elements(runs)


  ;; Get unique runs/reruns
  ten = ulong64(10)
  tid = ulong64(reruns)*ten^12 + ulong64(runs)*ten^15
  rmd = rem_dup(tid)
  Uruns = runs[rmd]
  Ureruns = reruns[rmd]

  ;; Will test whether between is faster than
  ;; specifying run and rerun

  where_statements = $
    "run = "+ntostr(Uruns)+" AND "+$
    "rerun = "+ntostr(Ureruns)


  struct = self->send_queries(table,where_statements, $
                              columns=columns, $
                              clauses=clauses, $
                              connect_info=connect_info, $
                              slow=slow, $
                              status=status)

  return,struct

END 


FUNCTION postgres_sdss::_get_reruns, runs, reruns=reruns, minrerun=minrerun

  nruns = n_elements(runs)
  nreruns = n_elements(reruns)
  IF nreruns EQ 0 THEN BEGIN 
      use_reruns = !sdss->rerun(runs, min=minrerun) 
      IF use_reruns[0] EQ -1 THEN BEGIN 
          message,'Halting'
      ENDIF 
      return, use_reruns
  ENDIF ELSE BEGIN 
      IF nreruns NE nruns THEN BEGIN 
          IF nreruns EQ 1 THEN BEGIN 
              return, replicate(reruns[0], nruns)
          ENDIF ELSE BEGIN 
              message,'reruns must be same size as runs or scalar'
          ENDELSE 
      ENDIF ELSE BEGIN 
          return, reruns
      ENDELSE 
  ENDELSE 

END 

FUNCTION postgres_sdss::readbyid, $
                      table, $
                      runs, camcols, fields, ids, $
                      reruns=reruns, minrerun=minrerun, $
                      columns=columns, $
                      clauses=clauses, $
                      connect_info=connect_info, $
                      slow=slow
                      

  np = n_params()
  IF np EQ 0 THEN BEGIN 
      print,'st = ps->readbyid(table, runs, camcols, fields, ids, '
      print,'                  reruns=, '
      print,'                  columns=, clauses=, '
      print,'                  connect_info='
      print,'                  /slow)'
      return,-1
  ENDIF 


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make calls depending on the entered id information
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF np EQ 1 THEN BEGIN 
      ;; The user wants to get entire runs
      struct = self->read_runs(table, runs, $
                               reruns=reruns, minrerun=minrerun, $
                               columns=columns, $
                               clauses=clauses, $
                               connect_info=connect_info, $
                               slow=slow)
      return,struct
  ENDIF 

  IF np EQ 2 THEN BEGIN 
      ;; The user wants to get entire runs,camcols
      struct = self->read_camcols(table, runs, camcols, $
                                  reruns=reruns, minrerun=minrerun, $
                                  columns=columns, $
                                  clauses=clauses, $
                                  connect_info=connect_info, $
                                  slow=slow)
      return,struct
  ENDIF 

  IF np EQ 3 THEN BEGIN 
      ;; The user wants individual runs,camcols,fields
      struct = self->read_fields(table, runs, camcols, fields, $
                                 reruns=reruns, minrerun=minrerun, $
                                 columns=columns, $
                                 clauses=clauses, $
                                 connect_info=connect_info, $
                                 slow=slow)
      return,struct
  ENDIF 

  IF np EQ 4 THEN BEGIN 
      ;; The user wants individual objects
      struct = self->read_ids(table, runs, camcols, fields, ids, $
                              reruns=reruns, minrerun=minrerun, $
                              columns=columns, $
                              clauses=clauses, $
                              connect_info=connect_info, $
                              slow=slow)
      return,struct
  ENDIF 

END 

PRO postgres_sdss__define

  struct = { postgres_sdss, $
             pgdummy: 0, $
             INHERITS postgres $
           }

END 
