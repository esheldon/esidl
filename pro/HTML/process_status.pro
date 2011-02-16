PRO process_status

  ;; make a web page with status of processing for each
  ;; stripe

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Fill in the comments for each stripe below
  ;; comments
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  comm_st = 37
  comments = 'Missing Chunk from end of run 1350'

  ;;;;;;;;;;;;;;;;;;;;;;
  ;; Stripe Index
  ;;;;;;;;;;;;;;;;;;;;;;

  read_stripe_specgal_inventory,invst

  ;;;;;;;;;;;;;;;;;;;;;;;
  ;; Output file
  ;;;;;;;;;;;;;;;;;;;;;;;

  outdir = '/net/cheops1/data0/esheldon/WWW/'
  outfile = outdir + 'completed_stripes.html'

  ;;;;;;;;;;;;;;;;;;;
  ;; get stripes
  ;;;;;;;;;;;;;;;;;;;

  run_status = sdss_runstatus() 
  w=where(run_status.stripe GT 0)

  run_status = run_status[w]

  rmd = rem_dup(run_status.stripe)
  stripes = run_status[rmd].stripe
  nst = n_elements(stripes)

  openw, lun, outfile, /get_lun

  printf, lun, '<html>'
  printf, lun, '<head>'
  printf, lun, '<title>Completed Stripes</title>'
  printf, lun, '</head>'
  printf, lun, '<body bgcolor="#ffffff" link="#0066ff" vlink="#009999" text="#000000">'
  printf, lun, '<h1>Stripe Information</h1>'
  printf, lun, '<ul>'
  printf, lun, '<li>The following table lists the processed runs for each stripe on the cheops cluster.'
  printf, lun, '</li>'
  printf, lun, '<li>Click on the stripe number to view the region of sky covered by the runs we have in the stripe (click on last column to see Steve Kent'+"'"+'s plots for all runs)'
  printf, lun, '<li>Click on the run number to view detailed info about the processing of that run/rerun'
  printf, lun, '</li>'
  printf, lun, '</li>'
  printf, lun, '<li>You can get this stripe info in IDL by typing sdss_goodstripe,stripe,runs,reruns.  This will print the runs/reruns for that stripe that are on the cheops cluster.'
  printf, lun, '</li>'
  printf, lun, '</ul>'
  printf, lun, 'Here is a polar <a href="./striperange/stripe_range.png">projection</a> of the survey area covered.'
  printf, lun, '<p><table border=1>'
  printf, lun, '<tr><th>Stripe</th><th>Strip-Run-Rerun</th><th>Photoz</th><th># Spectro Gal</th><th>Overall Comments</th><th>Steve Kent'+"'"+'s Plot</th></tr>'

  FOR i=0L, nst-1 DO BEGIN 

      ;; each element of the table
      stripe = stripes[i]
      ststr = ntostr(stripe)
      w = where(run_status.stripe EQ stripe, nrun)
      s = sort(run_status[w].strip)
      w=w[s]
      runs   = run_status[w].run
      reruns = run_status[w].rerun
      strips = run_status[w].strip

      ;; match to comments
      wcomm = where(comm_st EQ stripe, nw)
      IF nw EQ 0 THEN comm = '&nbsp' ELSE comm = comments[wcomm]

      ;; # of spectro galaxies
      winv = where(invst.stripe EQ stripe, ninv)
      IF ninv EQ 0 THEN ngalst = '0' ELSE ngalst = ntostr(invst[winv].ngal)

      col1 = '<a href="./striperange/stripe'+ststr+'_range.png">'+ststr+'</a>'
      col2 = ''
      col3 = ''
      FOR ir = 0L,nrun-1 DO BEGIN 
          rst = ntostr(runs[ir])
          rerst = ntostr(reruns[ir])
          stripst = strips[ir]
          col2 = col2 + '<a href="./process_status/'+rst+'/'+rerst+'/run'+rst+'-rerun'+rerst+'.html">'+stripst+'-'+rst+'-'+rerst+'</a><br>'

          ;; check for photoz
          si = sdss_flag_select(run_status.flags, 'runstatus', $
              {photoz_exist:'Y'}, input_index = w[ir])

          IF si[0] EQ -1 THEN hasphz = 'No' ELSE hasphz='Yes'
          IF ir LT nrun-1 THEN hasphz = hasphz + '<br>'
          col3 = col3 + hasphz
      ENDFOR
      col4 = ngalst
      col5 = comm
      col6 = $
        '<a href="http://www-sdss.fnal.gov:8000/skent/'+ststr+'-N.gif">'+ststr+'-N</a>&nbsp'+$
        '<a href="http://www-sdss.fnal.gov:8000/skent/'+ststr+'-S.gif">'+ststr+'-S</a>'
      printf,lun,$
             '<tr>'+$
             '<td nowrap align=center>'+col1+'</td>' + $
             '<td nowrap align=left>'+col2+'</td>' + $
             '<td nowrap align=left>'+col3+'</td>' + $
             '<td nowrap align=left>'+col4+'</td>' + $
             '<td nowrap align=left>'+col5+'</td>' + $
             '<td nowrap align=center>'+col6+'</td></tr>'

  ENDFOR 

  printf, lun, '</table>'
  printf, lun, '<hr>'
  printf, lun, '<b>Email: esheldon at cfcp.uchicago.edu</b>'
  printf, lun, 'Last modified: '+systime()
  printf, lun, '</body>'
  printf, lun, '</html>'

  free_lun, lun

END 
