PRO make_stripe_index, runs, stripes, strips

  ;; This makes a list or run/stripe/strip for all the
  ;; runs we have which have spectroscopy (in the directory
  ;; below)

  ;; This is useful because can be sent to programs
  ;; like run_combine_stripe_spec

  dir=sdssidl_config('SHAPECORR_DIR')+'spec_index/'
  outfile = dir+'stripe_index.dat'
  cd,dir

  files = findfile('run*')
  nf = n_elements(files)

  read_runstripe, rst
  read_runstripe, rst_local, /local

  wrst = where(rst.stripe GT 0 AND rst.stripe LE 86, nrst)
  rmd = rem_dup(rst[wrst].run)
  wrst = wrst[rmd]

  wrst_local = where(rst_local.stripe GT 0 AND rst_local.stripe LE 86, nrst_local)
  rmd = rem_dup(rst_local[wrst_local].run)
  wrst_local = wrst_local[rmd]

  runs = replicate(-1, nf)
  stripes = runs
  strips = replicate('?', nf)

  ;; get run numbers and their stripe info
  FOR i=0L, nf-1 DO BEGIN 

      f = files[i]
      tt = (str_sep(f,'_'))[0]
      tt2 = (str_sep(tt,'run'))[1]
      runs[i] = long(tt2)

      ;; match up
      w=where(rst[wrst].run EQ runs[i], nw)
      IF nw EQ 0 THEN BEGIN 
          print,'Checking the local fermi runs file: run '+ntostr(runs[i])
          w=where(rst_local[wrst_local].run EQ runs[i], nw)
          IF nw EQ 0 THEN BEGIN 

              print,'No match for run: '+ntostr(runs[i]),/inf
              tstripe = 0
              tstrip=' '
              read,tstripe, prompt='Enter stripe: '
              read,tstrip, prompt='Enter strip: '
              stripes[i] = tstripe
              strips[i] = tstrip
          ENDIF ELSE BEGIN 
              stripes[i] = rst_local[wrst_local[w]].stripe
              strips[i] = rst_local[wrst_local[w]].strip
          ENDELSE 
      ENDIF ELSE BEGIN 
          stripes[i] = rst[wrst[w]].stripe
          strips[i] = rst[wrst[w]].strip
      ENDELSE 

      print,runs[i],stripes[i],'  '+strips[i]

  ENDFOR 

  s = sort(runs)
  colprint, runs[s], stripes[s], '  '+strips[s], file=outfile

END 
