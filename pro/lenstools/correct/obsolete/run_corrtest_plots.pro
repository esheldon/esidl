PRO run_corrtest_plots, hirata=hirata

  ;; redo all them
  w=where(!run_status.adatc_photo_v ge 5.3, nw)
  w=where(!run_status.adatc_photo_v GE 5.3 AND !run_status.run EQ 756, nw)


  FOR i=0L, nw-1 DO BEGIN 
      run = !run_status[w[i]].run
      rerun = !run_status[w[i]].rerun

      FOR camcol=1,6 DO BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; regenerate files
          ;;;;;;;;;;;;;;;;;;;;;;;;;

          FOR clr=1,3 DO BEGIN 
              corrtest_plots, run, rerun, camcol, clr, hirata=hirata
          ENDFOR

      ENDFOR 

  ENDFOR 

END 
