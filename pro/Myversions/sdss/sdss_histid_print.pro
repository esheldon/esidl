PRO sdss_histid_print, baseStruct

  print,'Nobj:    '+ntostr(baseStruct.nobj)
  print,'Depth:   '+ntostr(basestruct.depth)
  print,'Nleaves: '+ntostr(baseStruct.nleaves)

  nruns = baseStruct.nruns
  pruns = baseStruct.runs
  FOR ri=0L, nruns-1 DO BEGIN

      run = (*pruns)[ri].run
      nreruns = (*pruns)[ri].nreruns
      IF nreruns GT 0 THEN BEGIN 

          preruns = (*pruns)[ri].reruns
          FOR rri=0L, nreruns-1 DO BEGIN 

              rerun = (*preruns)[rri].rerun

              ncamcols = (*preruns)[rri].ncamcols
              IF ncamcols GT 0 THEN BEGIN 
                  
                  pcamcols = (*preruns)[rri].camcols
                  FOR ci=0L, ncamcols-1 DO BEGIN 
                      camcol = (*pcamcols)[ci].camcol
                      
                      nfields = (*pcamcols)[ci].nfields

                      IF nfields GT 0 THEN BEGIN 

                          pfields = (*pcamcols)[ci].fields
                          FOR fi=0L, nfields-1 DO BEGIN 
                              field = (*pfields)[fi].field

                              nind = (*pfields)[fi].nind

                              print,run,rerun,camcol,field,nind
                          ENDFOR 

                      ENDIF ELSE BEGIN 
                          print,run,rerun,camcol
                      ENDELSE 
                  ENDFOR 

              ENDIF ELSE BEGIN 
                  print,run,rerun
              ENDELSE 

          ENDFOR 

      ENDIF ELSE BEGIN 
          print,run
      ENDELSE 
  ENDFOR 

END 
