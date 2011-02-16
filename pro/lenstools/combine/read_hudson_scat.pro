PRO read_hudson_scat, scat, rand

  stripes = [9,10,11,12,13,14,15]
  stripestr = stripearr2string(stripes)

  outdir = sdssidl_config('shapecorr_dir')+'combined/'
  srcfile = outdir + 'stripe'+stripestr+'_hudson.dat'

  instruct = create_struct('clambda', 0d, $
                           'ceta',    0d, $
                           'e1_recorr', 0.0, $
                           'e2_recorr', 0.0, $
                           'e1e1err',   0.0, $
                           'e1e2err',   0.0, $
                           'e2e2err',   0.0, $
                           'photoz_z',  0.0, $
                           'photoz_zerr', 0.0, $
                           'photoz_abscounts_u', 0.0, $
                           'photoz_abscounts_g', 0.0, $
                           'photoz_abscounts_r', 0.0, $
                           'photoz_abscounts_i', 0.0, $
                           'photoz_abscounts_z', 0.0, $
                           'upetro', 0.0, $
                           'gpetro', 0.0, $
                           'rpetro', 0.0, $
                           'ipetro', 0.0, $
                           'zpetro', 0.0, $
                           'rpetrorad', 0.0)

  print,'Counting lines in: ',srcfile
  nl = numlines(srcfile)
  scat = replicate(instruct, nl)

  print,'Reading file'
  openr, lun, srcfile, /get_lun
  tt = systime(1)
  readf, lun, scat
  free_lun, lun
  ptime,systime(1)-tt

  IF n_params() EQ 2 THEN BEGIN 

      print,'Reading rand file: ',rsrcfile
      rsrcfile = outdir + 'stripe'+stripestr+'_rand_hudson.dat'
      rinstruct =  create_struct('clambda', 0d, $
                                'ceta',    0d, $
                                'id', 0L)

      nl = numlines(rsrcfile)
      rand = replicate(rinstruct, nl)
      
      openr, lun, rinstruct, /get_lun
      readf, lun, rinstruct
      free_lun, lun

;      rdfloatstr, rsrcfile, rinstruct, rand, /double

  ENDIF 

END 
