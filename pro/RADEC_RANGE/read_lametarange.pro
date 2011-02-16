PRO read_lametarange, run, camcol, lametarange

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: read_lametarange, run, camcol, lametarange'
      return
  ENDIF 

  runstr = run2string(run)
  camstr = ntostr(camcol)

  file=sdssidl_config('radec_search_dir')+$
    'lameta-range-'+runstr+'-'+camstr+'.fit'

  lametarange = mrdfits(file,1)

END 
