;+
; NAME:
;  SDSS_MYSQL_READBYID
;
;
; PURPOSE:
;  Read objects from the database by id numbers.  The user can any number of
;  entire runs, or grab individual entries, or anything in between.  This
;  program optimizes the query for the most efficient search and retrieval of
;  rows. 
;
;
; CATEGORY:
;  SDSS routine.
;
;
; CALLING SEQUENCE:
;  sdss_mysql_readbyid, runs, camcols, fields, ids, 
;                       reruns=, /minrerun, 
;                       host=, user=, pwd=, database=, table=, 
;                       columns=, 
;                       /fast
;
;
; INPUTS:
;  runs: The runs to get.  This is the minimal number of input parameters. This
;        can be an array, in which case all unique runs in the array will be
;        read. 
;  camcols:  If runs and camcols are sent, all the unique run-camcol
;           combinations will be read.  The user can enter a single run with
;           multiple camcols, or a set of run-camcol pairs. 
;  fields: If runs,camcols,field sent, all unique sets are returned.  The user
;          can send a single pair of run,camcol and a set of fields, or the
;          full set of runs,camcols,fields can be sent.
;  ids:  If the user sends all id info, each array must be of the same length.
;
;
; OPTIONAL INPUTS:
;  reruns=: If not set, the default rerun is used, which is the latest unless
;           the /minrerun keyword is set, in which case the earliest rerun is
;           used. 
;  host=: The host machine.  By default this is retrieved from the sdssidl
;         configuration. 
;  user: The username. By default this is "sdss", for the anonymous, read-only
;        account.  
;  pwd: Password.  By default, it is assumed there is no password.
;  database=: which database to use. By default, it is the "sdss" database.
;  table=: Which table to use.  By default "tsObj"
;  columns=: Which columns to read for retrieved rows.  By default all are
;            read. (get a list of columns?)
;
;
; KEYWORD PARAMETERS:
;  /minrerun: If the user didn't enter the rerun info, use the
;             earliest. Default is to use the latest.
;  /fast: Use the "fast" method for retrieving rows.  
;
; OUTPUTS:
;  struct: A structure containing the results of the query.
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-

; Execute the input set of queries.

FUNCTION sdss_mysql_readbyid_queries, database, table, where_statements, $
                              host=host, user=user, pwd=pwd, $
                              columns=columns, $
                              fast=fast, $
                              nrows=nrows, status=status

  IF n_elements(columns) EQ 0 THEN BEGIN 
      colstr = '*'
  ENDIF ELSE BEGIN 
      colstr = strjoin(columns,', ')
  ENDELSE 
  
  nQuery = n_elements(where_statements)

  ;; Simple for a single query
  IF nQuery EQ 1 THEN BEGIN 

      query = "SELECT "+colstr+" from "+table+" WHERE "+where_statements[0]

      struct = mysql_query(database, query, $
                           host=host, user=user, pwd=pwd, $
                           status=status, nrows=nrows)
      return,struct

  ENDIF 

  ;; Multiple Queries sent
  IF keyword_set(fast) THEN BEGIN 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Fast: we gather all the objects in pointers and then
      ;; copy them into an output structure.  This is faster than the 
      ;; standard method below, but can use more memory for many separate
      ;; queries: a factor of two in the limit of an infinite number of
      ;; queries. 
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      ptrlist = ptrarr(nQuery)
      nrows = 0ULL
      
      FOR i=0L, nQuery-1 DO BEGIN 

          query = "SELECT "+colstr+" from "+table+$
            " WHERE "+where_statements[i]
          st = mysql_query(database, query, $
                           host=host, user=user, pwd=pwd, $
                           status=tstatus, nrows=tnrows)
          IF tnrows GT 0 THEN BEGIN 
              ptrlist[i] = ptr_new(st, /no_copy)
              nrows = nrows + tnrows
          ENDIF 

      ENDFOR 
      IF nrows GT 0 THEN BEGIN 
          struct = combine_ptrlist(ptrlist)
      ENDIF ELSE BEGIN 
          struct = -1
      ENDELSE 

  ENDIF ELSE BEGIN

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Standard: Count the queries results first by running a 
      ;; SELECT count(*) statement first, then create output and
      ;; copy in as we re-run queries.  For many large queries this 
      ;; could save lots of memory, but could also be a lot slower. 
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      numlist = ulon64arr(nQuery)
      nrows = 0ULL
      
      FOR i=0L, nQuery-1 DO BEGIN 

          query = "SELECT count(*) from "+table+$
            " WHERE "+where_statements[i]
          tnrows = mysql_query(database, query, $
                               host=host, user=user, pwd=pwd, $
                               status=tstatus, /noconvert) 

          IF strnumber(tnrows) THEN BEGIN 
              tnrows = ulong64(tnrows)
              nrows = nrows + tnrows
              numlist[i] = tnrows
          ENDIF 
          
      ENDFOR 

      ;; Any found?
      w=where(numlist GT 0, nw)
      IF nw NE 0 THEN BEGIN 

          ;; Now get a single structure
          query = "SELECT "+colstr+" from "+table+$
            " WHERE "+where_statements[w[0]]+" LIMIT 1"
          st =  mysql_query(database, query, $
                            host=host, user=user, pwd=pwd, $
                            status=tstatus) 

          ;; Now all
          beg = 0ULL
          struct = replicate(st, nrows)
          FOR i=0ULL, nQuery-1 DO BEGIN 

              IF numlist[i] NE 0 THEN BEGIN 

                  query = "SELECT "+colstr+" from "+table+$
                    " WHERE "+where_statements[i]
                  st = mysql_query(database, query, $
                                   host=host, user=user, pwd=pwd, $
                                   status=tstatus) 
                  
                  struct[beg:beg+numlist[i]-1] = temporary(st)
                  beg=beg+numlist[i]
                  
              ENDIF 

          ENDFOR 
          
      ENDIF ELSE BEGIN 
          ;; None found
          struct = -1
      ENDELSE 
      
  ENDELSE ;; Slower, memory efficient way

  return,struct

END 

;; The user has input the full id info for a list of objects

FUNCTION sdss_mysql_read_ids, database, table, runs, camcols, fields, ids, $
                              reruns=reruns, minrerun=minrerun, $
                              host=host, user=user, pwd=pwd, $
                              columns=columns, $
                              fast=fast, $
                              nrows=nrows, status=status

  nruns = n_elements(runs)
  nreruns = n_elements(reruns)
  ncamcols = n_elements(camcols)
  nfields = n_elements(fields)
  
  ;; Get the reruns if not sent
  IF nreruns EQ 0 THEN BEGIN 
      reruns=sdss_rerun(runs, min=minrerun) 
      nreruns = nruns

  ENDIF ELSE BEGIN 
      IF nreruns NE nruns THEN BEGIN 
          message,'reruns must be same size as runs',/inf
          return,-1
      ENDIF 
  ENDELSE 

  photoid = sdss_photoid(runs,reruns,camcols,fields,ids)

  photoid = photoid[rem_dup(photoid)]
  photoidstr = ntostr(photoid)
  photoidstr = '('+strjoin(photoidstr, ',')+')'

  where_statement = 'photoid IN '+photoidstr

  struct = sdss_mysql_readbyid_queries(database, table, where_statement, $
                                       host=host, user=user, pwd=pwd, $
                                       columns=columns, $
                                       fast=fast, $
                                       nrows=nrows, status=status)

  return,struct

END 

;; User has requested objects at the field level

FUNCTION sdss_mysql_read_fields, database, table, runs, camcols, fields, $
                                 reruns=reruns, minrerun=minrerun, $
                                 host=host, user=user, pwd=pwd, $
                                 columns=columns, $
                                 fast=fast, $
                                 nrows=nrows, status=status

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; For non-consecutive fields, /fast gives a 30% improvement for 100 fields
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  status = 1
  IF n_params() LT 2 THEN BEGIN 
      message,'write a syntax statement'
  ENDIF 

  nruns = n_elements(runs)
  nreruns = n_elements(reruns)
  ncamcols = n_elements(camcols)
  nfields = n_elements(fields)
  
  ;; Get the reruns if not sent
  IF nreruns EQ 0 THEN BEGIN 
      reruns=sdss_rerun(runs, min=minrerun) 
      nreruns = nruns

  ENDIF ELSE BEGIN 
      IF nreruns NE nruns THEN BEGIN 
          message,'reruns must be same size as runs',/inf
          return,-1
      ENDIF 
  ENDELSE 


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

      IF 1 THEN BEGIN 

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

          ;; Dumb way, slower unless all fields are non-consecutive
          where_statements = $
            "run = "+ntostr(runs)+" AND " + $
            "rerun = "+ntostr(reruns)+" AND "+$
            "camcol = "+ntostr(camcols)+" AND "+ $
            "field = "+ntostr(fields)

      ENDELSE 

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

      tp = 1
      CASE tp OF 
          0: BEGIN 

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
              ENDFOR            ;runs
              
              ;; Now combine the lists
              where_statements = combine_ptrlist(ptrlist)
              sdss_histid_destroy, idstruct
          END 
          1: BEGIN 
              ;; Get unique runs/reruns/camcols
              tt = replicate(1, nruns)
              id = sdss_photoid(runs,reruns,camcols,tt,tt)
          
              rmd = rem_dup(id)
              id = id[rmd]
              
              nid = n_elements(id)
              ptrlist = ptrarr(nid)
              ptrIndex = 0L
              
              ;; Histogram the runs
              Hrun = histogram(runs, rev=rrev)
              wHrun = where(Hrun NE 0, nwHrun)
              
              FOR ri=0L, nwHrun-1 DO BEGIN 
                  
                  ;; indices of this run in histogram
                  wHruni = wHrun[ri]
                  
                  ;; Indices of requests with this run
                  wrun = rrev[ rrev[wHruni]:rrev[wHruni+1] -1 ]
                  
                  ;; Histogram the reruns
                  Hrerun = histogram(reruns[wrun], rev=rrrev)
                  wHrerun = where(Hrerun NE 0, nwHrerun)
                  
                  FOR rri=0L, nwHrerun-1 DO BEGIN 
                      
                      ;; Indices of this rerun in rerun histogram
                      wHreruni = wHrerun[rri]
                      
                      ;; Indices of requests with this run-rerun
                      wrerun = rrrev[ rrrev[wHreruni]:rrrev[wHreruni+1]-1]
                      wrerun = wrun[wrerun]
                      
                      ;; Histogram the camcols
                      Hcamcol = histogram(camcols[wrerun], rev=crev)
                      wHcamcol = where(Hcamcol NE 0, nwHcamcol)
                      
                      FOR ci=0L, nwHcamcol-1 DO BEGIN 
                          
                          ;; Indices of this camcol in camcol histogram
                          wHcamcoli = wHcamcol[ci]
                          
                          ;; Indices of requests with this run-rerun-camcol
                          wcamcol = crev[ crev[wHcamcoli]:crev[wHcamcol+1] -1]
                          wcamcol = wrerun[wcamcol]
                          
                          ;; Now we are down to the different fields in this 
                          ;; camcol.  Same code as above, should make modular
                          
                          Ufields = fields[rem_dup(fields[wcamcol])]
                          nUfields = n_elements(Ufields)
                          
                          ;; Extract the consecutive ones
                          get_consecutive, uFields, first, last
                          nConsec = n_elements(first)
                          
                          ;; Generate where statements
                          twhere_statements = strarr(nConsec)
                          twhere_statements = $
                            "run = "+ntostr(runs[wcamcol])+" AND " + $
                            "rerun = "+ntostr(reruns[wcamcol])+" AND "+$
                            "camcol = "+ntostr(camcols[wcamcol])+" AND "+ $
                            "field BETWEEN "+$
                            ntostr(uFields[first])+" AND "+ntostr(uFields[last])
                          
                          ;; Now copy these into our pointerlist
                          ptrlist[ptrIndex] = ptr_new(twhere_statements, /no_copy)
                          ptrIndex = ptrIndex + 1
                          
                      ENDFOR ;; Camcols
                  ENDFOR ;; Reruns
              ENDFOR ;; Runs
              
              
              ;; Now combine the lists
              where_statements = combine_ptrlist(ptrlist)
          END 
          2: BEGIN 
              
              ;; Dumb way, doesn't take advantage of speedups
              print,'Doing straight'
              
              ;; Get unique runs/reruns/camcols/fields
              tt = replicate(1, nruns)
              id = sdss_photoid(runs,reruns,camcols,fields,tt)
              
              rmd = rem_dup(id)
              
              where_statements = $
                "run = "+ntostr(runs[rmd])+" AND "+$
                "rerun = "+ntostr(reruns[rmd])+" AND "+$
                "camcol = "+ntostr(camcols[rmd])+" AND "+$
                "field = "+ntostr(fields[rmd])
              where_statements = where_statements[ rem_dup(where_statements) ]

          END 
      ENDCASE 

  ENDELSE ;; multiple runs....

  struct = sdss_mysql_readbyid_queries(database, table, where_statements, $
                                       host=host, user=user, pwd=pwd, $
                                       columns=columns, $
                                       fast=fast, $
                                       nrows=nrows, status=status)

  return,struct

END 

;; The user wants entire camcols

FUNCTION sdss_mysql_read_camcols, database, table, runs, camcols, $
                                  reruns=reruns, minrerun=minrerun, $
                                  host=host, user=user, pwd=pwd, $
                                  columns=columns, $
                                  fast=fast, $
                                  nrows=nrows, status=status
  
  status = 1
  IF n_params() LT 2 THEN BEGIN 
      message,'write a syntax statement'
  ENDIF 

  nruns = n_elements(runs)
  nreruns = n_elements(reruns)
  ncamcols = n_elements(camcols)

  ;; Get the reruns if not sent
  IF nreruns EQ 0 THEN BEGIN 
      reruns=sdss_rerun(runs, min=minrerun) 
      nreruns = nruns
  ENDIF ELSE BEGIN 
      IF nreruns NE nruns THEN BEGIN 
          message,'reruns must be same size as runs',/inf
          return,-1
      ENDIF 
  ENDELSE 

  IF nruns EQ 1 THEN BEGIN 
      ;; Single run, perhaps multiple camcols
      where_statements = $
        "run = "+ntostr(runs)+" AND "+$
        "rerun = "+ntostr(reruns)+" AND "+$
        "camcol = "+ntostr(camcols)

  ENDIF ELSE BEGIN 

      ;; Multiple runs.  User must input exactly same number of camcols as
      ;; runs in this situation

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
  
  struct = sdss_mysql_readbyid_queries(database, table, where_statements, $
                                       host=host, user=user, pwd=pwd, $
                                       columns=columns, $
                                       fast=fast, $
                                       nrows=nrows, status=status)

  return,struct

END 

;; Entire runs are requested.

FUNCTION sdss_mysql_read_runs, database, table, runs, $
                               reruns=reruns, minrerun=minrerun, $
                               host=host, user=user, pwd=pwd, $
                               columns=columns, $
                               fast=fast, $
                               nrows=nrows, status=status
  status=1
  IF n_params() LT 1 THEN BEGIN 
      message,'write a syntax statement'
  ENDIF 

  nruns = n_elements(runs)
  nreruns = n_elements(reruns)

  ;; Get the reruns if not sent
  IF nreruns EQ 0 THEN BEGIN 
      reruns=sdss_rerun(runs, min=minrerun) 
      nreruns = nruns
  ENDIF ELSE BEGIN 
      IF nreruns NE nruns THEN BEGIN 
          message,'reruns must be same size as runs',/inf
          return,-1
      ENDIF 
  ENDELSE 


  ;; Get unique runs/reruns/camcols
  tt = replicate(1, nruns)
  id = sdss_photoid(runs,reruns,tt,tt,tt,tt)
  rmd = rem_dup(id)
  Uruns = runs[rmd]
  Ureruns = reruns[rmd]

  where_statements = $
    "run = "+ntostr(Uruns)+" AND "+$
    "rerun = "+ntostr(Ureruns)

  struct = sdss_mysql_readbyid_queries(database, table, where_statements, $
                                       host=host, user=user, pwd=pwd, $
                                       columns=columns, $
                                       fast=fast, $
                                       nrows=nrows, status=status)

  return,struct

END 


FUNCTION sdss_mysql_readbyid, runs, camcols, fields, ids, $
                              reruns=reruns, minrerun=minrerun, $
                              columns=columns, $
                              host=host, user=user, pwd=pwd, $
                              database=database, table=table, $
                              fast=fast

  np = n_params()

  IF np EQ 0 THEN BEGIN 
      print,'-Syntax: sdss_dbread, runs, camcols, fields, ids, '
      print,'            reruns=, /minrerun, '
      print,'            host=, user=, pwd=, database=, table=, '
      print,'            columns=, '
      print,'            /fast'
      return,-1
  ENDIF 

  ;; Read in the selected fields

  ;; Will change this to adatc or sdss or whatever.
  IF n_elements(database) EQ 0 THEN database = "sdss"
  IF n_elements(table) EQ 0 THEN table = "tsObj"

  IF np EQ 1 THEN BEGIN 
      ;; The user wants to get entire runs
      struct = sdss_mysql_read_runs(database, table, runs, $
                                    reruns=reruns, minrerun=minrerun, $
                                    host=host, user=user, pwd=pwd, $
                                    columns=columns, $
                                    fast=fast, $
                                    nrows=nrows, status=status)
      return,struct
  ENDIF 

  IF np EQ 2 THEN BEGIN 

      ;; The user wants entire camcols
      struct = sdss_mysql_read_camcols(database, table, runs, camcols, $
                                       reruns=reruns, minrerun=minrerun, $
                                       host=host, user=user, pwd=pwd, $
                                       columns=columns, $
                                       fast=fast, $
                                       nrows=nrows, status=status)

      return,struct

  ENDIF 

  IF np EQ 3 THEN BEGIN 
      ;; The user has specified fields as well

      struct = sdss_mysql_read_fields(database, table, runs, camcols, fields, $
                                      reruns=reruns, minrerun=minrerun, $
                                      host=host, user=user, pwd=pwd, $
                                      columns=columns, $
                                      fast=fast, $
                                      nrows=nrows, status=status)
      return,struct

  ENDIF 

  IF np EQ 4 THEN BEGIN 
      ;; id information is given: the user is requesting
      ;; individual objects

      struct = sdss_mysql_read_ids(database, table, $
                                   runs, camcols, fields, ids, $
                                   reruns=reruns, minrerun=minrerun, $
                                   host=host, user=user, pwd=pwd, $
                                   columns=columns, $
                                   fast=fast, $
                                   nrows=nrows, status=status)

      return,struct
  ENDIF 



END 
