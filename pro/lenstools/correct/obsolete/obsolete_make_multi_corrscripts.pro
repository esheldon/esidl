PRO print_shortstatus_check, lun

  printf, lun, 'if [ $? -ne 0 ]'
  printf, lun, 'then'
  printf, lun, '    echo ERROR: Did not finish'
  printf, lun, 'fi'
  printf, lun


END 

PRO make_multi_corrscripts, stripe, color_index, outfront, stripeind=stripeind

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: make_multi_corrscripts, stripe, color_index [, outfront, stripeind=stripeind]'
      print,'stripeind for subset of runs/reruns'
      print,'outfront to override default corr_stripe09 type front'
      return
  ENDIF 

  get_good_lenstripe, stripe, runs, reruns, psFieldreruns=psFieldreruns

  nrun=n_elements(runs)
  FOR i=0L, nrun-1 DO BEGIN 
      make_corrscripts, runs[i], reruns[i], color_index, $
                        psFieldrerun=psFieldreruns[i], $
                        script_dir=sdir, outf_dir=odir

      add_arrval, temporary(sdir), script_dirs
      add_arrval, temporary(odir), outf_dirs
  ENDFOR 

  print
  IF n_elements(outfront) EQ 0 THEN BEGIN 
      outfront = 'corr_stripe'+stripe2string(stripe)
  ENDIF
  
  outdir = !sdss_shapecorr_dir + 'scripts/corrstripe/'
  outfile = outdir + outfront+'-col1-2-3.sh'

  print,'Making script: ',outfile
  openw, lun, outfile, /get_lun

  printf, lun, '#!/bin/sh'
  FOR i=0L, nrun-1 DO BEGIN 
      rstr=ntostr(runs[i])
      printf, lun
      printf, lun, 'echo Processing run '+rstr
      printf, lun, $
              script_dirs[i]+'correct'+rstr+'-col1-2-3.sh'+$
              ' 1> '+ outf_dirs[i]+'correct'+rstr+'-col1-2-3.out' + $
              ' 2> '+ outf_dirs[i]+'correct'+rstr+'-col1-2-3.err'
      print_shortstatus_check, lun
  ENDFOR 
  free_lun, lun
  spawn,'chmod 755 '+outfile

  print
  outfile = outdir + outfront+'-col4-5-6.sh'
  print,'Making script: ',outfile
  openw, lun, outfile, /get_lun
  printf, lun, '#!/bin/sh'
  FOR i=0L, nrun-1 DO BEGIN 
      rstr=ntostr(runs[i])
      printf, lun
      printf, lun, 'echo Processing run '+rstr
      printf, lun, $
              script_dirs[i]+'correct'+rstr+'-col4-5-6.sh'+$
              ' 1> '+ outf_dirs[i]+'correct'+rstr+'-col4-5-6.out' + $
              ' 2> '+ outf_dirs[i]+'correct'+rstr+'-col4-5-6.err'
      print_shortstatus_check, lun
  ENDFOR 
  free_lun, lun
  spawn,'chmod 755 '+outfile

END 
