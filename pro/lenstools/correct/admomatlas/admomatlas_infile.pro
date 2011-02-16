PRO admomatlas_infile, run, rerun, camcol, fname_in, fname_out

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: admomatlas_infile, run, rerun, camcol, fname_in, fname_out'
      return
  ENDIF 

  fname_in = 'admomin-'+run2string(run) + '-' + $
                        ntostr(camcol) + '-' + $
                        ntostr(rerun) + '.dat'

  fname_out = repstr(fname_in, 'admomin', 'admomout')

   return

END 
