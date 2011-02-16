PRO run_corrtest_plots, corrdirs=corrdirs,$
                        remove=remove, no_overwrite=no_overwrite

  ;; redo all the corrtest stuff for every run on current current cheops 
  ;; machine
  ;; can override this with corrdirs=corrdirs
  ;;
  ;; set /remove to remove all corshape* files before running corrtest_plots
  ;; set /overwrite to overwrite ps files

  ;; camcols/colors to do

  colors=[1,2,3] & nclr=n_elements(colors)
  camcols=[1,2,3,4,5,6] & ncamcol=n_elements(camcols)

  IF n_elements(corrdirs) EQ 0 THEN BEGIN 
      corrdirs = ["/data0/corrected.local/", $
                 "/data1/corrected.local/"]
  ENDIF 

  print,'Corrdirs: '
  forprint,corrdirs

  print
  ncorr=n_elements(corrdirs)

  FOR icorr=0L, ncorr-1 DO BEGIN 
      corrdir = corrdirs[icorr]

      ;; try to change to that dir
      !error_state.name = ""
      cd,corrdir

      ;; try to change dirs
      IF !error_state.name NE "IDL_M_CNTCNGDIR" THEN BEGIN 

          print
          print,"Processing corrdir: ",corrdir
          print

          ;; get all the corrdirs
          spawn,'ls -d corr[1-9]*',dirlist
          IF dirlist[0] EQ '' THEN BEGIN
              ;; not need to crash here, its ok if nothing here
              print,"no corrdirs in this directory: "+corrdir
              return
          ENDIF 
          
          print
          print,"Processing these runs: "
          forprint,dirlist
          print
          
          ;; Ok, now corrtest stuff for each 
          nrun = n_elements(dirlist)
          FOR i=0L, nrun-1 DO BEGIN 
              
              ;; find reruns
              run = long( repstr( dirlist[i], "corr", "") )
              
              w=where(!run_status.run EQ run, nrerun)
              
              ;; loop over the reruns
              FOR j=0L, nrerun-1 DO BEGIN 
                  
                  rerun = !run_status[w[j]].rerun
                  FOR cam=0L, ncamcol-1 DO BEGIN 
                      
                      camcol = camcols[cam]
                      
                      ;; remove files if requested
                      IF keyword_set(remove) THEN BEGIN 
                          command = "rm corr"+ntostr(run)+"/"+ntostr(rerun)+$
                            "/objcs/"+ntostr(camcol)+"/corshape*"
                          print,command
                          spawn,command
                      ENDIF 
                      FOR clr=0L, nclr-1 DO BEGIN 
                          color=colors[clr]
                          
                          
                          corrtest_plots, run, rerun, camcol, color,$
                                          no_overwrite=no_overwrite
                          
                      ENDFOR ;; clr
                  ENDFOR ;; camcol
              ENDFOR ;; rerun
          ENDFOR ;; run
      ENDIF ELSE BEGIN 
          print,"Can't chdir dir to "+corrdir
      ENDELSE 

  ENDFOR ;; corrdirs loop
END 
