PRO shapelet_decomp_scripts, set, typestr

  dir = '/net/cheops2/home/esheldon/idlscripts/shapelet/'

  setstr = ntostr(set)
  CASE set OF
      1: BEGIN 
          ;; Already did 9
;          stripes = [9,10,11,12,13,14,15,16]
          stripes = [10,11,12,13,14,15,16]
      END 
      2: BEGIN 
          stripes=[26,27,28,29,30,31,32,33]
      END 
      3: BEGIN 
          stripes=[34,35,36,37,42,43,44,76,82,86]
      END 
      ELSE: message,'Unknown set #: '+setstr
  ENDCASE 

  file = dir + 'run_shapelet_decomp_'+typestr+'_'+setstr+'.sh'
  openw, lun, file, /get_lun

  printf, lun, '!/bin/sh'
  printf, lun

  nstripe = n_elements(stripes)
  FOR i=0L, nstripe-1 DO BEGIN 

      sstr = ntostr(stripes[i])

      printf, lun, 'nice -19 idl<<EOF'
      printf, lun, '  stripe='+sstr
      printf, lun, '  type="'+typestr+'"'
      printf, lun, '  nmax=15'
      printf, lun, '  clr=2'
      printf, lun, '  run_shapelet_decomp, stripe, type, nmax, clr, $'
      printf, lun, '          /doplot, /zbuff, /png, status=status'
      printf, lun, '  if status ne 0 then exit, status=45'
      printf, lun, 'EOF'
      printf, lun, 'status=$?'
      printf, lun, 'if [ $status -ne 0 ]'
      printf, lun, 'then'
      printf, lun, '    echo Error in `basename $0` for stripe '+sstr
      printf, lun, '    err="Shapelet Error for stripe '+sstr+'"'
      printf, lun, '    echo "$err" | mail esheldon@cfcp.uchicago.edu -s "$err"'
      printf, lun, 'fi'
      printf, lun
  ENDFOR 

  printf,lun,'dt=`date`'
  printf,lun,'message="Finished shapelets set '+setstr+'"'
  printf,lun,'echo "$message:  $dt" | mail esheldon@cfcp.uchicago.edu -s "$message"'

  free_lun, lun

  spawn,['chmod','755',file],/noshell

END 
