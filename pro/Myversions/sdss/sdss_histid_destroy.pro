PRO sdss_histid_destroy, baseStruct

  ;; We allow indices to be stored at all levels.  Check base structure.  
  IF ptr_valid(baseStruct.indices) THEN BEGIN 
      ptr_free, baseStruct.indices
  ENDIF 

  pruns = baseStruct.runs
  IF ptr_valid(pruns) THEN BEGIN 
      
      ;; Now loop over runs
      nruns = baseStruct.nruns
      FOR ri=0L, nruns-1 DO BEGIN 

          ;; Check for indices on basic run structure
          IF ptr_valid( (*pruns)[ri].indices ) THEN BEGIN 
              ptr_free, (*pruns)[ri].indices 
          ENDIF 

          ;; Reruns for this run?
          preruns = (*pruns)[ri].reruns
          IF ptr_valid(preruns) THEN BEGIN 

              ;; Now loop over the reruns
              nreruns = (*pruns)[ri].nreruns
              FOR rri=0L, nreruns-1 DO BEGIN 

                  ;; Check for indices on the base rerun structure
                  IF ptr_valid( (*preruns)[rri].indices ) THEN BEGIN 
                      ptr_free, (*preruns)[rri].indices 
                  ENDIF 
                  
                  ;; Camcols for this rerun?
                  pcamcols = (*preruns)[rri].camcols
                  IF ptr_valid(pcamcols) THEN BEGIN 

                      ;; Now loop over the camcols
                      ncamcols = (*preruns)[rri].ncamcols
                      FOR ci=0L, ncamcols-1 DO BEGIN 

                          ;; Check for indices on the base camcol structure
                          IF ptr_valid( (*pcamcols)[ci].indices ) THEN BEGIN 
                              ptr_free, (*pcamcols)[ci].indices 
                          ENDIF 

                          ;; Fields for this camcol?
                          pfields = (*pcamcols)[ci].fields
                          IF ptr_valid(pfields) THEN BEGIN 
                              
                              nfields = (*pcamcols)[ci].nfields
                              FOR fi=0L, nfields-1 DO BEGIN 

                                  ;; Check for indices for this field
                                  pind = (*pfields)[fi].indices
                                  IF ptr_valid(pind) THEN BEGIN 
                                      ptr_free, pind
                                  ENDIF 
                              ENDFOR 
                              ptr_free, pfields
                          ENDIF ;; There are fields for this camcol

                      ENDFOR 
                      ptr_free,pcamcols
                  ENDIF ;; There are camcols for this reruns
                  
              ENDFOR 
              ptr_free, preruns
          ENDIF ;; There are reruns for this run

      ENDFOR 
      ptr_free, pruns

  ENDIF ;; There are individual runs

END 
