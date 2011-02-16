PRO get_meanscinv, stripe, clr, scinv_struct, $
                   use_lambda=use_lambda, hirata=hirata,$
                   lrg_sources=lrg_sources, rlrg_sources=rlrg_sources

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: get_meanscinv, stripe, clr, scinv_struct, /use_lambda, /hirata, lrg_sources=lrg_sources'
      return
  ENDIF 

  meanscinv_names, stripe, clr, fitfile, $
    use_lambda=use_lambda, hirata=hirata, $
    lrg_sources=lrg_sources, rlrg_sources=rlrg_sources

  IF NOT fexist(fitFile) THEN BEGIN 
      message,'File does not exist: '+fitFile
  ENDIF 

  print
  print,'Reading mean scinv struct: ',fitfile
 
  scinv_struct = mrdfits(fitfile, 1)

END 
