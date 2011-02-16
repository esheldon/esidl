PRO run_tests, gal=gal

  ;; compare directly
  outdir = '/net/cheops2/home/esheldon/tmp/'

  colmin = 1
  colmax = 6
  clrmin = 1
  clrmax = 3

  ;;GOTO, jump1

  runs = [745,756]
  nrun =  n_elements(runs)
  FOR irun=0L, nrun-1 DO BEGIN 
      run=runs[irun]
      FOR camcol = colmin,colmax DO BEGIN 
          FOR clr = clrmin,clrmax DO BEGIN 
              
              ;; comparison w/old
              IF keyword_set(gal) THEN addstr='gal_' ELSE addstr=''
              file = outdir + 'compare_'+addstr+'newold'+ntostr(run)+'-'+ntostr(camcol)+'-'+!colors[clr]+'.ps'
              begplot, name=file, /color
              test_vsold_samefile, struct, run, camcol, clr, gal=gal
              endplot

              ;; psf correction
              ;file = outdir + 'psfcorr_newold'+ntostr(run)+'-'+ntostr(camcol)+'-'+!colors[clr]+'.ps'
              ;begplot, name=file, /color
              ;test_psfcorr, struct, run, camcol, clr
              ;endplot


          ENDFOR
          delvarx, struct
      ENDFOR
  ENDFOR 
  return
jump1:
  ;; errors
  FOR camcol = colmin,colmax DO BEGIN 
      FOR clr = clrmin,clrmax DO BEGIN 

          file = outdir + 'testerror_745_756-'+ntostr(camcol)+'-'+!colors[clr]+'.ps'
          begplot, name=file, /color
          test_error, prun1, prun2, camcol, clr
          endplot

      ENDFOR
      delvarx, prun1, prun2
  ENDFOR


END 
  
