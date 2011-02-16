PRO print_shortstatus_check, lun

  printf, lun, 'if [ $? -ne 0 ]'
  printf, lun, 'then'
  printf, lun, '    echo ERROR: Did not finish'
  printf, lun, 'fi'
  printf, lun


END 

PRO make_multi_corrscripts_byrun, runs, reruns, outfront, maxmag, disk=disk, addphotoz=addphotoz, addbayes=addbayes, nonewcrun=nonewcrun

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: make_multi_corrscripts_byrun, runs, reruns, outfront, maxmag, [disk=disk, addphotoz=addphotoz, addbayes=addbayes]'
      return
  ENDIF 

  nrun=n_elements(runs)
  ndisk = n_elements(disk)
  
  FOR i=0L, nrun-1 DO BEGIN 
      
      IF ndisk NE 0 THEN BEGIN 
          IF ndisk EQ 1 THEN sdisk=disk[0] $
          ELSE sdisk = disk[i]
      ENDIF 
      
      make_corrscripts, runs[i], reruns[i], maxmag, $
        scriptDir=sdir, outFileDir=odir, disk=sdisk,$
        addphotoz=addphotoz, addbayes=addbayes, $
        nonewcrun=nonewcrun
      
      add_arrval, temporary(sdir), scriptDirs
      add_arrval, temporary(odir), outFileDirs
  ENDFOR 

  spawn,'addrunlinks corrected'

  outdir = sdssidl_config('shapecorr_dir') + 'scripts/multirun/'
  outfile = outdir + outfront+'-col1-2-3.sh'

  print,'Making script: ',outfile
  openw, lun, outfile, /get_lun


  printf, lun, '#!/bin/sh'
  FOR i=0L, nrun-1 DO BEGIN 
      rstr  = ntostr(runs[i])
      rrstr = ntostr(reruns[i])
      printf, lun
      printf, lun, 'echo Processing run '+rstr+' rerun '+rrstr
      printf, lun, $
              scriptDirs[i]+'correct'+rstr+'-col1-2-3.sh'+$
              ' 1> '+ outFileDirs[i]+'correct'+rstr+'-col1-2-3.out' + $
              ' 2> '+ outFileDirs[i]+'correct'+rstr+'-col1-2-3.err'
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
      rstr  = ntostr(runs[i])
      rrstr = ntostr(reruns[i])
      printf, lun
      printf, lun, 'echo Processing run '+rstr+' rerun '+rrstr
      printf, lun, $
              scriptDirs[i]+'correct'+rstr+'-col4-5-6.sh'+$
              ' 1> '+ outFileDirs[i]+'correct'+rstr+'-col4-5-6.out' + $
              ' 2> '+ outFileDirs[i]+'correct'+rstr+'-col4-5-6.err'
      print_shortstatus_check, lun
  ENDFOR 
  free_lun, lun
  spawn,'chmod 755 '+outfile

END 
