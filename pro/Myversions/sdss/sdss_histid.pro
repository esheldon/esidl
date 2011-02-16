FUNCTION sdss_histid_histogram, data, rev=rev

  IF n_elements(data) EQ 1 THEN BEGIN 
      hist = histogram([data], rev=rev)
  ENDIF ELSE BEGIN 
      hist = histogram(data, rev=rev)
  ENDELSE 

  return, hist
END 

FUNCTION sdss_histid_checkpars, np, runs, reruns, camcols, fields

  nruns = n_elements(runs)
  nreruns = n_elements(reruns)
  ncamcols = n_elements(camcols)
  nfields = n_elements(fields)

  CASE np OF
      1: return,1
      2: BEGIN 
          IF nreruns NE nruns THEN BEGIN 
              message,'runs,reruns must be same length',/inf
              return,0 
          ENDIF ELSE return,1
      END 
      3: BEGIN 
          IF nreruns NE nruns OR ncamcols NE nruns THEN BEGIN
              message,'runs,reruns,camcols must be same length',/inf
              return,0 
          ENDIF ELSE return,1
      END 
      4: BEGIN 
          IF (nreruns NE nruns OR $
              ncamcols NE nruns OR $
              nfields NE nruns) THEN BEGIN 
              message,'runs,reruns,camcols,fields must be same length',/inf
              return,0
          ENDIF ELSE return,1
      END  
      ELSE:
  ENDCASE 

END 

FUNCTION sdss_histid, runs, reruns, camcols, fields, $
                      status=status, check=check

  status = 1
  np = n_params()

  IF np LT 1 OR np GT 4 THEN BEGIN 
      print,'-Syntax: idstruct = sdss_histid(runs, [reruns, camcols, fields, status=])'
      return,-1
  ENDIF 

  IF NOT sdss_histid_checkpars(np, runs, reruns, camcols, fields) THEN BEGIN 
      return,-1
  ENDIF 

  ;; Node structures. 
  ;; Depth     meaning
  ;;  1    only runs indexed
  ;;  2    runs/reruns
  ;;  3    runs/reruns/camcols
  ;;  4    runs/reruns/camcols/fields

  ;; indices will remain null except at the leaf of the tree
  baseStruct   = {depth:0, nleaves:0L, nobj:0L, $
                  nruns:0, runs: ptr_new(), nind:0L, indices:ptr_new()}

  runStruct    = {run:0L, $
                  nreruns:0,  reruns:  ptr_new(), nind:0L, indices:ptr_new()}
  rerunStruct  = {rerun:0,$
                  ncamcols:0, camcols: ptr_new(), nind:0L, indices:ptr_new()}
  camcolStruct = {camcol:0,  $
                  nfields:0,  fields:  ptr_new(), nind:0L, indices:ptr_new()}
  fieldStruct  = {field:0, nind:0L, indices:ptr_new()}

  ;; Histogram the runs
  run_h = sdss_histid_histogram(runs, rev=rrev)
  wrun_h = where(run_h NE 0, nwrun_h)
  
  ;; Fill in the main run struct
  baseStruct.depth = np
  baseStruct.nruns = nwrun_h
  baseStruct.runs = ptr_new(replicate(runStruct, nwrun_h))
  pruns = baseStruct.runs

  FOR ri=0L, nwrun_h-1 DO BEGIN 
      
      ;; indices of this run in histogram
      wrun_h_i = wrun_h[ri]
      
      ;; Indices of requests with this run
      wrun = rrev[ rrev[wrun_h_i]:rrev[wrun_h_i+1] -1 ]

      ;; Copy in run
      (*pruns)[ri].run = runs[wrun[0]]

      ;; Are there reruns?
      IF np GT 1 THEN BEGIN 
          ;; Histogram the reruns
          rerun_h = sdss_histid_histogram(reruns[wrun], rev=rrrev)
          wrerun_h = where(rerun_h NE 0, nwrerun_h)
          
          ;; set rerun structures, and get a pointer to them
          (*pruns)[ri].nreruns = nwrerun_h
          (*pruns)[ri].reruns = $
            ptr_new(replicate(rerunStruct, nwrerun_h))
          
          preruns = (*pruns)[ri].reruns
          FOR rri=0L, nwrerun_h-1 DO BEGIN 

              ;; Indices of this rerun in rerun histogram
              wrerun_h_i = wrerun_h[rri]
              
              ;; Indices of requests with this run-rerun
              wrerun = rrrev[ rrrev[wrerun_h_i]:rrrev[wrerun_h_i+1]-1]
              wrerun = wrun[wrerun]

              ;; Copy in rerun
              (*preruns)[rri].rerun = reruns[wrerun[0]]

              ;; Are there camcols?
              IF np GT 2 THEN BEGIN               
                  ;; Histogram the camcols
                  hamcol_h = sdss_histid_histogram(camcols[wrerun], rev=crev)
                  wCamcol_h = where(hamcol_h NE 0, nwcamcol_h)
                  
                  ;; Set camcols structures, and get a pointer to them
                  (*preruns)[rri].ncamcols = nwcamcol_h
                  (*preruns)[rri].camcols = $
                    ptr_new(replicate(camcolStruct, nwcamcol_h))

                  pcamcols = (*preruns)[rri].camcols                  
                  FOR ci=0L, nwcamcol_h-1 DO BEGIN 
                      
                      ;; Indices of this camcol in camcol histogram
                      wcamcol_h_i = wcamcol_h[ci]
                      
                      ;; Indices of requests with this run-rerun-camcol
                      wcamcol = crev[ crev[wcamcol_h_i]:crev[wcamcol_h_i+1] -1]
                      wcamcol = wrerun[wcamcol]

                      ;; Copy in camcol
                      (*pcamcols)[ci].camcol = camcols[wcamcol[0]]

                      ;; Are there fields?
                      IF np GT 3 THEN BEGIN 
                      
                          ;; Histogram the fields
                          field_h = $
                            sdss_histid_histogram(fields[wcamcol], rev=frev)
                          wfield_h = where(field_h NE 0, nwfield_h)
                          
                          ;; Set fields structures, and get a pointer to them
                          (*pcamcols)[ci].nfields = nwfield_h
                          (*pcamcols)[ci].fields = $
                            ptr_new(replicate(fieldStruct,nwfield_h))

                          pfields = (*pcamcols)[ci].fields
                          FOR fi=0L, nwfield_h-1 DO BEGIN 
                              
                              ;; Indices of this field in field histogram
                              wfield_h_i = wfield_h[fi]
                              
                              ;; indices of requests with this 
                              ;; run-rerun-camcol-field
                              wfield = $
                                frev[ frev[wfield_h_i]:frev[wfield_h_i+1] -1 ]
                              wfield = wcamcol[wfield]

                              nind = n_elements(wfield)
                              (*pfields)[fi].field = fields[wfield[0]]
                              (*pfields)[fi].indices = ptr_new(wfield,/no_copy)
                              (*pfields)[fi].nind = nind

                              baseStruct.nleaves = baseStruct.nleaves+1
                              baseStruct.nobj = baseStruct.nobj + nind
                              baseStruct.nleaves = baseStruct.nleaves + 1
                          ENDFOR ;; Fields
                      ENDIF ELSE BEGIN ;; np < 3
                          nind = n_elements(wcamcol)
                          (*pcamcols)[ci].indices = ptr_new(wcamcol, /no_copy)
                          (*pcamcols)[ci].nind = nind
                          baseStruct.nleaves = baseStruct.nleaves+1
                          baseStruct.nobj = baseStruct.nobj + nind
                          baseStruct.nleaves = baseStruct.nleaves + 1
                      ENDELSE 
                  ENDFOR ;; Camcols
              ENDIF ELSE BEGIN ;; np < 2
                  nind = n_elements(wrerun)
                  (*preruns)[rri].indices = ptr_new(wrerun, /no_copy)
                  (*preruns)[rri].nind = nind

                  baseStruct.nobj = baseStruct.nobj + nind
                  baseStruct.nleaves = baseStruct.nleaves + 1
              ENDELSE 
          ENDFOR ;; Reruns
      ENDIF ELSE BEGIN ;; np < 1
          nind = n_elements(wrun)
          (*pruns)[ri].indices = ptr_new(wrun, /no_copy)
          (*pruns)[ri].nind = nind

          baseStruct.nobj = baseStruct.nobj + nind
          baseStruct.nleaves = baseStruct.nleaves + 1
      ENDELSE 
  ENDFOR ;; Runs

  IF keyword_set(check) THEN BEGIN 

      ;; Check total number of indices, as well as match the ids up to the
      ;; indexed objects

      ntotal = 0ULL

      print,'Checking'
      nn = n_elements(wstruct)
      FOR i=0LL, nn-1 DO BEGIN 
          
          ind = *wstruct[i].indices
          
          ntotal = ntotal + n_elements(ind)

          w=where(runs[ind] NE wstruct[i].run, nw)
          IF nw NE 0 THEN message,'runs'

          IF np GT 1 THEN BEGIN 
              w=where(reruns[ind] NE wstruct[i].rerun, nw)
              IF nw NE 0 THEN message,'reruns'
              
              IF np GT 2 THEN BEGIN 
                  w=where(camcols[ind] NE wstruct[i].camcol, nw)
                  IF nw NE 0 THEN message,'camcols'
                  
                  IF np GT 3 THEN BEGIN 
                      w=where(fields[ind] NE wstruct[i].field, nw)
                      IF nw NE 0 THEN message,'fields'
                  ENDIF ;; fields
              ENDIF ;; camcols
          ENDIF ;; reruns
          ;; runs

      ENDFOR 

      IF ntotal NE n_elements(runs) THEN BEGIN 
          message,'Total number of indices not equal input'
      ENDIF 

  ENDIF 

  status = 0
  return,baseStruct


END 
