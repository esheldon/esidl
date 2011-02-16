PRO run_check_admom, runs, reruns, bad, badfiles

  nrerun = n_elements(reruns)
  bad = bytarr(nrerun, 6)
  badfiles = ptrarr(nrerun, 6)
  
  totstr = ntostr(nrerun)
  FOR i=0L, nrerun-1 DO BEGIN 

      run = runs[i]
      rerun = reruns[i]

      print
      print,'Checking: '+ntostr(run)+' '+ntostr(rerun)+$
            ' ('+ntostr(i+1)+'/'+totstr+')'
      FOR camcol=1,6 DO BEGIN 
          print,'Camcol: '+ntostr(camcol)+'  ',format='(a,$)'
          check_admom, run, rerun, camcol, tbad
          IF n_elements(tbad) NE 0 THEN BEGIN 
              print,'Found Bad'
              bad[i, camcol-1] = 1b
              badfiles[i, camcol-1] = ptr_new(tbad, /no_copy)
          ENDIF ELSE BEGIN 
              print,'Ok'
          ENDELSE 

      ENDFOR 

      print,'Hit a key'
      key=get_kbrd(1)

  ENDFOR 


END 
