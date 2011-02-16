PRO split_sxwtheta_byrerun, file

  ;; split sx output into files by run/rerun
  ;; then use combine_sxwtheta to combine into stripes

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: split_sxwtheta_byrerun, file'
      return
  ENDIF 

  print
  print,'Input file  ',file
  print
  tf=str_sep(file, '.fit')

  tf=tf[0]

  print,'Front ends begin with ',tf

  str=mrdfits(file,1)

  runs=str[rem_dup(str.run)].run
  nrun=n_elements(runs)

  outst = create_struct('run', 0L, $
                        'gr', 0.0,$
                        'petrocounts',fltarr(5),$
                        'lambda',0d,$
                        'eta',0d)

  gr = ( (str.counts_model[1] - str.reddening[1]) - $
         (str.counts_model[2] - str.reddening[2]) )
  petrocounts = str.petrocounts - str.reddening

  FOR irun=0L, nrun-1 DO BEGIN 

      run=runs[irun]

      wrun=where(str.run EQ run,nrun)
      rdup=rem_dup(str[wrun].rerun)
      reruns = str[wrun[rdup]].rerun
      nrerun = n_elements(reruns)

      FOR irer=0L, nrerun-1 DO BEGIN

          rerun=reruns[irer]

          wrer=where(str[wrun].rerun EQ rerun, nrer)
          print,'Run: ',ntostr(run),'  Rerun: ',ntostr(rerun),' Nobj: ',ntostr(nrer)

          tmp = replicate(outst, nrer)
          
          tmp.run = run
          tmp.gr = gr[wrun[wrer]]
          tmp.petrocounts = petrocounts[*,wrun[wrer]]
          
          eq2survey, str[wrun[wrer]].ra, str[wrun[wrer]].dec,lam, eta
          tmp.lambda = lam
          tmp.eta = eta

          newfile = tf+'-'+run2string(run)+'-'+ntostr(rerun)+'.fit'
          print,'Output file: ',newfile

          hdr0=['END']
          sxaddpar, hdr0, 'run', run
          sxaddpar, hdr0, 'rerun', rerun

          mwrfits2, tmp, newfile, /create, /destroy, hdr0=hdr0

          delvarx, lam, eta

      ENDFOR 

  ENDFOR 

  delvarx,str

END 
